import os
import subprocess
import time
import random
import glob
import re
from multiprocessing import Pool
import os.path

DEMAND_FRACTION = 0.05
TRY_COUNT = 1000
AGGREGATE_COUNT = 100000
OPT = "B4"

topologies = ["Arpanet196912.graph","Renam.graph","Mren.graph","Telecomserbia.graph","Basnet.graph","Nordu1989.graph","Padi.graph","Epoch.graph","Layer42.graph","Napnet.graph","Sanren.graph","Gblnet.graph","Getnet.graph","HiberniaIreland.graph","Netrail.graph","Ai3.graph","Cesnet1993.graph","Dataxchange.graph","JanetExternal.graph","Arpanet19706.graph","Nordu2005.graph","Nsfcnet.graph","Itnet.graph","Singaren.graph","Heanet.graph","Cesnet1997.graph","Cesnet1999.graph","Kreonet.graph","Abilene.graph","Nordu1997.graph","HiberniaCanada.graph","Nsfnet.graph","Ilan.graph","TLex.graph","Eenet.graph","HiberniaUk.graph","Nordu2010.graph","Grena.graph","Navigata.graph","Spiralight.graph","Sprint.graph","Compuserve.graph","Gridnet.graph","KentmanJul2005.graph","Claranet.graph","Eunetworks.graph","Garr199901.graph","Jgn2Plus.graph","Marwan.graph","Rhnet.graph","Sago.graph","Peer1.graph","Nextgen.graph","HostwayInternational.graph","Fatman.graph","Savvis.graph","Twaren.graph","Restena.graph","Arpanet19719.graph","HiberniaNireland.graph","Uninet.graph","Atmnet.graph","BsonetEurope.graph","Ibm.graph","Aarnet.graph","Ans.graph","Cesnet2001.graph","GtsRomania.graph","Istar.graph","Noel.graph","Garr200109.graph","Garr200404.graph","KentmanApr2007.graph","Renater1999.graph","VisionNet.graph","Harnet.graph","York.graph","Azrena.graph","Oxford.graph","Pacificwave.graph","Uran.graph","Amres.graph","Garr199904.graph","Garr199905.graph","Marnet.graph","Psinet.graph","Packetexchange.graph","Garr200112.graph","Bandcon.graph","Easynet.graph","Fccn.graph","Vinaren.graph","Renater2001.graph","Goodnet.graph","EliBackbone.graph","HiberniaUs.graph","Arpanet19723.graph","Globalcenter.graph","BtAsiaPac.graph","KentmanFeb2008.graph"]

exec_location = "../build/stability_eval"
topologies_location = "/Users/nik/Downloads/zoo_data/data/2016TopologyZooUCL_inverseCapacity"

def RunCmd(cmd):
#    print 'Will run', cmd
    subprocess.call(cmd, shell=True)
    #print 'Done'

def ProcessTopology(topology):
    top_location = '{}/{}'.format(topologies_location, topology)
    tm_string = '{}/{}.*.demands'.format(topologies_location, topology.split('.graph')[0])

    cmd = '{} --topology_file {} --matrices="{}" --demand_fraction {} --try_count {} --aggregate_count {} --opt {}'.format(exec_location, top_location, tm_string, DEMAND_FRACTION, TRY_COUNT, AGGREGATE_COUNT, OPT)
    RunCmd(cmd)

for top in topologies:
    ProcessTopology(top)
