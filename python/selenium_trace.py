from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time
import subprocess
import os
import signal

PAGES = ['http://www.google.co.uk', 'http://www.youtube.com',
         'http://www.google.com', 'http://www.facebook.com',
         'http://www.reddit.com', 'http://www.amazon.co.uk',
         'http://www.bbc.co.uk', 'http://www.wikipedia.org',
         'http://www.ebay.co.uk', 'http://www.twitter.com',
         'http://www.ladbible.com', 'http://www.live.com',
         'http://www.yahoo.com', 'http://www.instagram.com',
         'http://www.linkedin.com', 'http://www.livejasmin.com',
         'http://www.netflix.com', 'http://www.theguardian.com',
         'http://www.imgur.com', 'http://www.diply.com',
         'http://www.dailymail.co.uk', 'http://www.paypal.com',
         'http://www.vk.com', 'http://www.imdb.com',
         'http://www.twitch.tv', 'http://www.tumblr.com',
         'http://www.givemesport.com', 'http://www.gov.uk',
         'http://www.wikia.com', 'http://www.sportbible.com',
         'http://www.office.com', 'http://www.rightmove.co.uk',
         'http://www.booking.com', 'http://www.gumtree.com',
         'http://www.telegraph.co.uk', 'http://www.tripadvisor.co.uk',
         'http://www.msn.com', 'http://www.bing.com',
         'http://www.wordpress.com', 'http://www.microsoft.com',
         'http://www.gov.uk', 'http://www.lloydsbank.co.uk',
         'http://www.pirateproxy.cc', 'http://www.stackoverflow.com']
IFACE = 'eno1'

def RunWget(output):
    pcap_cmd = 'tcpdump -n -i {} -w {}'.format(IFACE, output)
    process = subprocess.Popen(
        pcap_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        shell=True, preexec_fn=os.setsid) 
    return process

def DoPages(driver, wait_time=2):
    for page in PAGES:
        driver.get(page)
        time.sleep(wait_time)

def AddAdblockToProfile(profile):
    profile.add_extension("adblock.xpi")
    profile.set_preference("extensions.adblockplus.suppress_first_run_page", "true")

def DisableCacheForProfile(profile):
    profile.set_preference("browser.cache.disk.enable", False)
    profile.set_preference("browser.cache.memory.enable", False)
    profile.set_preference("browser.cache.offline.enable", False)
    profile.set_preference("network.http.use-cache", False)

for use_adblock in [False]:
    # Will first do a run without the cache.
    profile = webdriver.FirefoxProfile()
    if use_adblock:
        AddAdblockToProfile(profile)
    DisableCacheForProfile(profile)

    driver = webdriver.Firefox(firefox_profile=profile)
    p = RunWget('adblock_{}_nocache'.format(use_adblock))
    time.sleep(2)

    DoPages(driver)
    os.killpg(p.pid, signal.SIGTERM)
    time.sleep(2)
    driver.close()
    
    # Now need to get a new profile, with the cache enabled, rerun the
    # pages to warm up the cache and then do another run.
    profile = webdriver.FirefoxProfile()
    if use_adblock:
        AddAdblockToProfile(profile)

    driver = webdriver.Firefox(firefox_profile=profile)
    DoPages(driver, wait_time=0)
    
    p = RunWget('adblock_{}_cache'.format(use_adblock))
    time.sleep(2)

    DoPages(driver)
    os.killpg(p.pid, signal.SIGTERM)
    time.sleep(2)
    driver.close()
