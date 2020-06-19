# -*- coding: utf-8 -*-
""" 
Python Script
Created on  Tuesday April 2019 09:04:53 
@author:  domy7 

[desc]
Description of the plugin Here
Write here any thing... 
[/desc]

ARGUMENTS:
----------
<inp> 
    _input :[required] - [type = int] - [default = None] 
    Descripe your input here 
        * bullet point.
        * bullet point
</inp>
<inp>
    Other inputs go here ...
</inp>

RETURN:
----------
    <out>
        output_ : indicate your output description here. \n refers to a new line.
    </out>

"""
import Grasshopper as gh
output_=gh.addPoint(0.,4.,2.)
a=gh.addPoint(-5.,10.,0.)

p1=gh.Point(0.,0.,0.)
p2=gh.Point(5.,5.,5.)

b=gh.addLine(p1,p2)