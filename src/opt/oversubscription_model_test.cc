#include "oversubscription_model.h"

#include "test_common.h"
#include "gtest/gtest.h"

namespace fubar {
namespace test {

using namespace ncode::net;

class ModelTest : public TestBase {
 protected:
  void AddAggregate(uint64_t cookie, const std::string& src,
                    const std::string& dst, Bandwidth total_volume,
                    size_t num_flows) {
    AggregateInput input({cookie, graph_storage_.NodeFromStringOrDie(src),
                          graph_storage_.NodeFromStringOrDie(dst)},
                         num_flows, total_volume);
    AggregateOutput output(input);
    aggregate_map_.emplace(cookie, output);
  }

  void AddPath(uint64_t cookie, const std::string& path_string,
               double fraction) {
    AggregateOutput& aggregate_output =
        ncode::FindOrDie(aggregate_map_, cookie);
    const GraphPath* path =
        graph_storage_.PathFromStringOrDie(path_string, cookie);

    PathOutput path_output(fraction, path);
    aggregate_output.AddToPaths(path_output);
  }

  void CheckPath(uint64_t cookie, const std::string& path_string,
                 Bandwidth per_flow_bps, Bandwidth bps_limit,
                 Bandwidth backup_bps) {
    const AggregateOutput& aggregate_output =
        ::ncode::FindOrDie(output_map_, cookie);
    const GraphPath* path =
        graph_storage_.PathFromStringOrDie(path_string, cookie);
    bool path_found = false;
    for (const PathOutput* path_output : aggregate_output.PathsSorted()) {
      if (path_output->path() == path) {
        ASSERT_FALSE(path_found);
        ASSERT_DOUBLE_EQ(path_output->per_flow_bps(), per_flow_bps.bps());
        ASSERT_DOUBLE_EQ(path_output->bps_limit(), bps_limit.bps());
        ASSERT_DOUBLE_EQ(path_output->backup_limit(), backup_bps.bps());
        path_found = true;
      }
    }
    ASSERT_TRUE(path_found);
  }

  size_t BackupPathCount(uint64_t cookie) {
    const AggregateOutput& aggregate_output =
        ::ncode::FindOrDie(output_map_, cookie);
    size_t i = 0;
    for (const PathOutput* path_output : aggregate_output.PathsSorted()) {
      if (path_output->backup_limit() > 0) {
        ++i;
      }
    }

    return i;
  }

  void CheckOversubscribedLink(const std::string& src, const std::string& dst,
                               double level) {
    const GraphStorage* graph_storage = path_cache_.graph_storage();
    bool link_found = false;
    for (const auto& link_index_and_info : oversubscribed_links_) {
      GraphLinkIndex link_index = link_index_and_info.first;
      const OutputLinkInfo& link_info = *link_index_and_info.second;

      const GraphLink* link_ptr = graph_storage->GetLink(link_index);
      if (link_ptr->src_node()->id() == src &&
          link_ptr->dst_node()->id() == dst) {
        ASSERT_FALSE(link_found);
        ASSERT_DOUBLE_EQ(link_info.subscription, level);
        link_found = true;
      }
    }
    ASSERT_TRUE(link_found);
  }

  void RunModel(double capacity_multiplier) {
    OverSubModel model(capacity_multiplier, aggregate_map_, &path_cache_);
    output_map_ = model.GetModifiedOutputMap();

    for (const auto& link_index_and_info : model.GetLinksInfo()) {
      GraphLinkIndex link_index = link_index_and_info.first;
      const OutputLinkInfo& link_info = *link_index_and_info.second;

      if (link_info.subscription > 1.0001) {
        oversubscribed_links_[link_index] = link_info;
      }
    }
  }

  std::map<uint64_t, AggregateOutput> aggregate_map_;
  std::map<uint64_t, AggregateOutput> output_map_;
  GraphLinkMap<OutputLinkInfo> oversubscribed_links_;
};

TEST_F(ModelTest, Empty) {
  OverSubModel model(1.0, aggregate_map_, &path_cache_);
  ASSERT_TRUE(model.GetModifiedOutputMap().empty());
  ASSERT_EQ(0ul, model.GetLinksInfo().Count());
}

TEST_F(ModelTest, SinglePathFit) {
  AddAggregate(100, "A", "C", kDefaultLinkspeed, 1000);
  AddPath(100, "[A->B, B->C]", 1.0);

  RunModel(1.0);

  ASSERT_EQ(1ul, output_map_.size());
  CheckPath(100, "[A->B, B->C]", kDefaultLinkspeed / 1000, kDefaultLinkspeed,
            Bandwidth::Zero());
  CheckPath(100, "[A->D, D->C]", Bandwidth::Zero(), Bandwidth::Zero(),
            kDefaultLinkspeed);
  ASSERT_EQ(0ul, oversubscribed_links_.Count());
}

TEST_F(ModelTest, SinglePathNoFit) {
  AddAggregate(100, "A", "C", kDefaultLinkspeed * 2, 1000);
  AddPath(100, "[A->B, B->C]", 1.0);

  RunModel(1.0);

  ASSERT_EQ(1ul, output_map_.size());
  CheckPath(100, "[A->B, B->C]", kDefaultLinkspeed / 1000, kDefaultLinkspeed,
            Bandwidth::Zero());
  CheckPath(100, "[A->D, D->C]", Bandwidth::Zero(), Bandwidth::Zero(),
            kDefaultLinkspeed);
  ASSERT_EQ(2ul, oversubscribed_links_.Count());
  CheckOversubscribedLink("A", "B", 2.0);
  CheckOversubscribedLink("B", "C", 2.0);
}

TEST_F(ModelTest, TwoPathsFit) {
  AddAggregate(100, "A", "C", kDefaultLinkspeed / 2.0, 1000);
  AddAggregate(101, "A", "D", kDefaultLinkspeed / 2.0, 1000);
  AddPath(100, "[A->B, B->C]", 1.0);
  AddPath(101, "[A->B, B->C, C->D]", 1.0);

  RunModel(1.0);

  ASSERT_EQ(2ul, output_map_.size());
  CheckPath(100, "[A->B, B->C]", kDefaultLinkspeed / 2000,
            kDefaultLinkspeed / 2, Bandwidth::Zero());
  CheckPath(101, "[A->B, B->C, C->D]", kDefaultLinkspeed / 2000,
            kDefaultLinkspeed / 2, Bandwidth::Zero());

  ASSERT_EQ(1ul, BackupPathCount(100));
  CheckPath(100, "[A->D, D->C]", Bandwidth::Zero(), Bandwidth::Zero(),
            kDefaultLinkspeed);
  ASSERT_EQ(1ul, BackupPathCount(101));
  CheckPath(101, "[A->D]", Bandwidth::Zero(), Bandwidth::Zero(),
            kDefaultLinkspeed);
  ASSERT_EQ(0ul, oversubscribed_links_.Count());
}

TEST_F(ModelTest, TwoPahsNoFitDifferentFlowCount) {
  AddAggregate(101, "A", "C", kDefaultLinkspeed / 2, 500);
  AddAggregate(100, "A", "D", kDefaultLinkspeed, 1000);
  AddPath(100, "[A->B, B->C, C->D]", 1.0);
  AddPath(101, "[A->B, B->C]", 1.0);

  RunModel(1.0);

  CheckPath(100, "[A->B, B->C, C->D]", kDefaultLinkspeed / 1000.0 / 1.5,
            kDefaultLinkspeed / 1.5, Bandwidth::Zero());
  CheckPath(101, "[A->B, B->C]", kDefaultLinkspeed / 2.0 / 500.0 / 1.5,
            kDefaultLinkspeed / 2.0 / 1.5, Bandwidth::Zero());
}

TEST_F(ModelTest, TwoPahsNoFitDifferentVolume) {
  AddAggregate(100, "A", "C", kDefaultLinkspeed * 10, 1000);
  AddAggregate(101, "A", "D", kDefaultLinkspeed * 2, 1000);
  AddPath(100, "[A->B, B->C]", 1.0);
  AddPath(101, "[A->B, B->C, C->D]", 1.0);

  RunModel(1.0);

  CheckPath(100, "[A->B, B->C]", (10 * kDefaultLinkspeed) / 1000.0 / 12,
            (10 * kDefaultLinkspeed) / 12.0, Bandwidth::Zero());
  CheckPath(101, "[A->B, B->C, C->D]", (2 * kDefaultLinkspeed) / 1000.0 / 12,
            (2 * kDefaultLinkspeed) / 12.0, Bandwidth::Zero());
}

TEST_F(ModelTest, TwoPathsFitNoBackup) {
  AddAggregate(100, "A", "C", kDefaultLinkspeed, 1000);
  AddAggregate(101, "A", "D", kDefaultLinkspeed, 1000);
  AddPath(100, "[A->D]", 1.0);
  AddPath(101, "[A->B, B->C]", 1.0);

  RunModel(1.0);

  ASSERT_EQ(2ul, output_map_.size());
  CheckPath(101, "[A->B, B->C]", kDefaultLinkspeed / 1000.0, kDefaultLinkspeed,
            Bandwidth::Zero());
  CheckPath(100, "[A->D]", kDefaultLinkspeed / 1000.0, kDefaultLinkspeed,
            Bandwidth::Zero());

  // There should be no backup paths, as the 2 paths have taken up all the
  // network.
  ASSERT_EQ(0ul, BackupPathCount(100));
  ASSERT_EQ(0ul, BackupPathCount(101));
  ASSERT_EQ(0ul, oversubscribed_links_.Count());
}

}  // namespace
}  // namespace fubarizer
