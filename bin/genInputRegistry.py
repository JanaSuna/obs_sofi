#!/usr/bin/env python2
# 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011, 2012, 2013 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import argparse
import glob
import os
import re
try:
    import sqlite3
except ImportError:
    # try external pysqlite package; deprecated
    import sqlite as sqlite3
import sys
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage

parser = argparse.ArgumentParser()
parser.add_argument("--create", default=False, action="store_true", help="Create new registry (clobber old)?")
parser.add_argument("--root", default=".", help="Root directory")
args = parser.parse_args()

root = args.root
calibs = glob.glob(os.path.join(root, "*","*.fits.gz"))
fields = glob.glob(os.path.join(root, "*","*","*.fits.gz"))
files = calibs+fields

registryName = os.path.join(root,"registry.sqlite3")
if os.path.exists(registryName) and args.create:
    os.unlink(registryName)

makeTables = not os.path.exists(registryName)
conn = sqlite3.connect(registryName)
if makeTables:
    cmd = "create table raw (id integer primary key autoincrement"
    cmd += ", expNum int, filename text, pointing text, filter text, dateObs text, mjd double, ra text, dec text, expTime double)"
    
    #mjd,pointing,dateObs,id,filter,filename,ra,expNum,dec,expTime
    
    conn.execute(cmd)
    conn.commit()
    cmd = "create table raw_visit (id integer primary key autoincrement"
    cmd += ", expNum int, filename text, pointing text, filter text, dateObs text, mjd double, ra text, dec text, expTime double)"
    conn.execute(cmd)
    conn.commit()
    cmd = "create table raw_flat (id integer primary key autoincrement"
    cmd += ", expNum int, filename text, pointing text, filter text, dateObs text, mjd double, expTime double)"
    conn.execute(cmd)
    cmd = "create table raw_dark (id integer primary key autoincrement"
    cmd += ", expNum int, filename text, pointing text, filter text, dateObs text, mjd double, expTime double)"
    conn.execute(cmd)
    conn.commit()

for fits in files:
    matches = re.search(r'(FIELDS|FLAT|DARK|STANDARD)/(((\w{9})/(\w{3})_(\w{3})_(\d{2}))|(FLAT_(\w{5}))|(STD_(\d{4})_(\w{5}))|D_\d{1}|D_\d{2}|D_1.2)_(\d{3}).fits', fits)
    #File names examples: FIELDS: FIELDS/05feb_F02_S22_10_021.fits
    #FIELDS/%(dateObs)s_F02_%(pointing)s_%(expTime)d_%(expNum)03d.fits
    #Flats: FLAT/FLAT_06Feb_005.fits
    #Standards: STANDARD/STD_9104_05feb_053.fits
    if not matches:
        print >>sys.stderr, "Warning: skipping unrecognized filename:", fits
        continue

    type = fits[len(root):len(root)+4]

    if (type=='FIEL'):
        dateObs = fits[len(root)+11:len(root)+11+5]
        pointing = fits[len(root)+21:len(root)+21+3]
    if (type=='FLAT'):
        dateObs = fits[len(root)+10:len(root)+10+5]
        print dateObs
    if (type=='DARK'):
        dateObs = 'Unknown'

    expNum = fits[len(fits)-8-3:len(fits)-8]


    # Extract information from header
    im = afwImage.ExposureF(fits)
    h = im.getMetadata()
    mjd = float(h.get('MJD-OBS'))
    filename = h.get('ARCFILE')
    objtype = h.get('OBJECT')
    filter = 'Ks'
    pointing = ''
    expTime = float(h.get('ESO DET DIT'))

    try:
        if type=='FIEL':
            ra = h.get('RA')
            dec= h.get('DEC')
            conn.execute("INSERT INTO raw VALUES (NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                     (expNum, filename, pointing, filter, dateObs, mjd, ra, dec, expTime))

        #conn.execute("INSERT INTO raw_visit VALUES (NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?)",                          (expNum, dateObs, mjd, expTime, ra, dec, filename, objtype, filter))
        if type=='FLAT':
            conn.execute("INSERT INTO raw_flat VALUES (NULL, ?, ?, ?, ?, ?, ?, ?)",
                     (expNum, filename, pointing, filter, dateObs, mjd, expTime))
        if type=='DARK':
            conn.execute("INSERT INTO raw_dark VALUES (NULL, ?, ?, ?, ?, ?, ?, ?)",
                     (expNum, filename, pointing, filter, dateObs, mjd, expTime))
    
    except Exception as e:
        print ("skipping botched %s: %s" % (fits, e))
        continue

conn.commit()
conn.close()
