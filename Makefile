#!/usr/bin/env make -f
#
# $Id$

VERSION=`svn info |grep Revision|sed -e 's/Revision: //g'`
DATE=`svn info |grep 'Last Changed Date: '|sed -e 's/Last Changed Date: //g'|cut -f1 -d' '`

TARGET=/Users/gaser/spm/spm8/toolbox/TOM

STARGET=141.35.200.101:/Applications/xampp/htdocs/

FILES=TOM.man tbx_cfg_tom.m cg_tom.m spm_TOM.m Contents.m Install.txt Changes

ZIPFILE=TOM8_$(VERSION).zip

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
	-@echo '% __________________________________________________________________________' > Contents.m
	-@echo '% Template-O-Matic Toolbox' >> Contents.m
	-@echo '% Version ' ${VERSION} ' (TOM8) ' ${DATE} >> Contents.m
	-@cat Contents_info.txt >> Contents.m

zip: update
	-@echo zip
	-@test ! -d TOM || rm -r TOM
	-@cp -rp ${TARGET} .
	-@zip ${ZIPFILE} -rm TOM

scp: zip
	-@echo scp
	-@scp -pr ${ZIPFILE} ${STARGET}
