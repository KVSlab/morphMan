from os import path, listdir
from IPython import embed
import importlib

files = [f[:-3] for f in listdir(".") if f[-3:] == ".py"]

functocheck = "get_teams"

"""
for f in files:
    testmod = importlib.import_module(f)
    func = dir(testmod)
    if functocheck in func:
        embed()
        print "ITS HERE"
        print 1/0
"""
    
from modulefinder import ModuleFinder

finder = ModuleFinder()
finder.run_script('move_siphon.py')

print 'Loaded modules:'
for name, mod in finder.modules.iteritems():
    print '%s: ' % name,
    print ','.join(mod.globalnames.keys()[:3])

print '-'*50
print '-'*50
print '-'*50
print '-'*50
print '-'*50
print '-'*50
print '-'*50
print '-'*50
print '-'*50
print '-'*50
print '-'*50
print '-'*50
print 'Modules not imported:'
print '\n'.join(finder.badmodules.iterkeys())
