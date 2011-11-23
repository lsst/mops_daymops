#!/usr/bin/env python


import sys

if __name__=="__main__":
     total = 0
     tr = 0 
     fa = 0
     
     lines = sys.stdin.readlines()

     for line in lines:
          items = line.split()
          if line[0] != "!":
               total += int(items[0])
               tr += int(items[1])
               fa += int(items[2])
     print "Total: ", total, " true Tracks: ", tr, " false tracks: ", fa
     if total > 0:
          print " %True: ", tr*100./total
               
