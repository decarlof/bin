#!/usr/bin/env python3

"""
Launch a script if specified files change.

"""

import os
import os.path
from pyinotify import WatchManager, IN_DELETE, IN_CREATE, IN_CLOSE_WRITE, IN_MOVED_TO, ProcessEvent, Notifier
import re
import subprocess
import sys
import argparse
import fnmatch

__author__ = "Alexander Bernauer (alex@copton.net) https://github.com/copton/react"
__copyright__ = "[GPL 2.0](http://www.gnu.org/licenses/gpl-2.0.html)"
__docformat__ = 'restructuredtext en'
__all__ = ['main']

class PatternAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, fnmatch.translate(values))

parser = argparse.ArgumentParser(description='Launch a script if specified files change.')
parser.add_argument('directory', help='the directory which is recursively monitored')

group = parser.add_mutually_exclusive_group()
group.add_argument('-r', '--regex', required=False, default=".*", help='files only trigger the reaction if their name matches this regular expression')
group.add_argument('-p', '--pattern', required=False, dest="regex", action=PatternAction, help='files only trigger the reaction if their name matches this shell pattern')

parser.add_argument("script", help="the script that is executed upon reaction")

class Options:
    __slots__=["directory", "regex", "script"]

options = Options()
args = parser.parse_args(namespace=options)

class Reload (Exception):
    pass

class Process(ProcessEvent):
    def __init__(self,  options):
        self.regex = re.compile(options.regex)
        self.script = options.script

    def process_IN_CREATE(self, event):
        target = os.path.join(event.path, event.name)
        if os.path.isdir(target):
            raise Reload()

    def process_IN_DELETE(self, event):
        raise Reload()

    # used for catching mv hidden to normal file
    def process_IN_MOVED_TO(self, event):
        target = os.path.join(event.path, event.name)
        if self.regex.match(target):
            args = self.script.replace('$f', target).split()
            #os.system("clear")
            sys.stdout.write("executing script: " + " ".join(args) + "\n")
            subprocess.call(args)
            sys.stdout.write("------------------------\n")

    def process_IN_CLOSE_WRITE(self, event):
        target = os.path.join(event.path, event.name)
        if self.regex.match(target):
            # Ignore all hidden '.*' files 
            if len(event.name) > 0:
                if event.name[0] == '.':
                    return
            args = self.script.replace('$f', target).split()
            #os.system("clear")
            sys.stdout.write("executing script: " + " ".join(args) + "\n")
            subprocess.call(args)
            sys.stdout.write("------------------------\n")

while True:
    wm = WatchManager()
    process = Process(options)
    notifier = Notifier(wm, process)
    mask = IN_DELETE | IN_CREATE | IN_CLOSE_WRITE | IN_MOVED_TO 
    wdd = wm.add_watch(options.directory, mask, rec=True)
    try:
        while True:
            notifier.process_events()
            if notifier.check_events():
                notifier.read_events()
    except Reload:
        pass
    except KeyboardInterrupt:
        notifier.stop()
        break
