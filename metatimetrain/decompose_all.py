#!/usr/bin/env python
# coding: utf-8
#[] wrapper to decompose iteratively all matrices.

from metatimetrain import decompose
def main( inputlist, thread):
    decompose.main('decompose 1')
    decompose.main('decompose 2')
    print(thread)
