"""Use the selected Python's approved image/numeric libraries before optional extras."""
import sys,runpy,argparse
import PIL.Image
import numpy
p=argparse.ArgumentParser();p.add_argument('--extra-packages');p.add_argument('script');p.add_argument('arguments',nargs=argparse.REMAINDER);a=p.parse_args()
if a.extra_packages:sys.path.insert(0,a.extra_packages)
sys.argv=[a.script]+a.arguments
runpy.run_path(a.script,run_name='__main__')
