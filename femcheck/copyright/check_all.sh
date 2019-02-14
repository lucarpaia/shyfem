#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# checks all folders for copyright
#
#---------------------------------------------------

basedir=$HOME/shyfem
copydir=$HOME/shyfem/femcheck/copyright
check=$copydir/check.log
check_orig=$copydir/check_orig.log
aux=$copydir/aux.log

#---------------------------------------------------

CheckDir()
{
  echo "========================================"
  echo "./$1"
  echo "========================================"
  $copydir/check_copyright_iter.sh $basedir/$1
}

CleanLog()
{
  $copydir/check_all.pl $1 $check > $aux
  mv $aux $check
}

ElabLog()
{
  echo "+++++++++++++++++++++++++++++++++++++++"
  cat $check | sort | uniq \
	| grep -v -- '-----' \
	| grep -v '=====' \
	| grep -v '^\.' \
	| grep -v '\.ps$' \
	| grep -v '\.eps$' \
	| grep -v '\.pdf$' \
	| grep -v '\.bat$' \
	| grep -v '\.grd$' 
  echo "+++++++++++++++++++++++++++++++++++++++"
}

#---------------------------------------------------

[ -f $check ] && rm $check

if [ 1 -eq 1 ]; then
CheckDir examples		| tee -a $check
#CheckDir fem3d			| tee -a $check
#CheckDir fembin			| tee -a $check
#CheckDir femregress		| tee -a $check
fi

if [ 1 -eq 0 ]; then
CheckDir femadj			| tee -a $check
CheckDir femanim		| tee -a $check
CheckDir femcheck		| tee -a $check
CheckDir femdoc			| tee -a $check
CheckDir femdummy		| tee -a $check
CheckDir femlib			| tee -a $check
CheckDir femplot		| tee -a $check
CheckDir femspline		| tee -a $check
CheckDir femutil		| tee -a $check
CheckDir grid			| tee -a $check
CheckDir hcbs			| tee -a $check
CheckDir mesh			| tee -a $check
CheckDir post			| tee -a $check
fi

cp $check $check_orig

CleanLog ./color
CleanLog ./perl/modules
CleanLog ./perl/codepage
CleanLog ./perl/GD
CleanLog ./python

ElabLog

#---------------------------------------------------

