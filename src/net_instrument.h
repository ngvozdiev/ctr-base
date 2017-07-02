#ifndef NET_INSTRUMENT
#define NET_INSTRUMENT

namespace ctr {

// Instruments the network by periodically recording per-queue and per-link
// stats to metrics.
class NetInstrument : public nc::EventConsumer {
 public:
  void HandleEvent() override;

 private:
  void Record();

  std::vector<const nc::htsim::Queue*> queues_;
  std::vector<const nc::htsim::Pipe*> pipes_;
};

}  // namespace ctr

#endif
