#!/usr/bin/env make -f
#
# $Id$

VERSION=`svn info |grep Revision|sed -e 's/Revision: //g'`
DATE=`svn info |grep 'Last Changed Date: '|sed -e 's/Last Changed Date: //g'|cut -f1 -d' '`

TARGET=/Users/gaser/spm/spm8/toolbox/TOM8

STARGET=141.35.200.101:/Applications/xampp/htdocs/TOM8

FILES=TOM8.man tbx_cfg_tom8.m cg_tom8.m cg_tom8_update.m spm_TOM8.m Contents.m Install.txt Readme-NIH-Data.txt Changes

ZIPFILE=TOM8_r$(VERSION).zip

install:
	-@echo install
	-@test ! -d ${TARGET} || rm -r ${TARGET}
	-@mkdir ${TARGET}
	-@cp ${FILES} ${TARGET}

help:
	-@echo Available commands:
	-@echo install zip scp update

update:
	-@svn update
	-@echo '% Template-O-Matic Toolbox' > Contents.m
	-@echo '% Version ' ${VERSION} ' (TOM8) ' ${DATE} >> Contents.m
	-@cat Contents_info.txt >> Contents.m

zip: update
	-@echo zip
	-@test ! -d TOM8 || rm -r TOM8
	-@cp -rp ${TARGET} .
	-@zip ${ZIPFILE} -rm TOM8

scp: zip
	-@echo scp to http://dbm.neuro.uni-jena.de/TOM8/${ZIPFILE}
	-@scp -pr ${ZIPFILE} ${STARGET}
