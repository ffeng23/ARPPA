# ARPPA
an R package for Protoarray Analysis
#README -- this file is not bundled together with package to distribute.
#This file is to describe miscellaneous things/considerations about building
#,updating and maintaining the package. it is more about programming 
#and housekeeping. 
#
#To read things more about theory and how the package works themselves,
#please go to Note.txt under ./dev folder

1.directory structure of the package
	./git -- git utility for storing necessary database, files, etc
	./data ./R ./man -- r package standard directory folders. ./R is the only
			mandatory folder that is created upon initialization by either
			skeleton() or setup() or create(). ./man is created only after
			calling document(). ./data could either created manually or
			by calling use_data() or use_data_raw(). Check the documentations
			for creating/adding/documenting data in 
			feng/docs/R_material/Advanced_R/"Data ¡¤ R packages.pdf"
			or go to url: http://rpkgs.had.co.nz/data.html
	./inst -- is also one the standard directory folders that will be used 
			by r package utility. The nice thing about this folder is that 
			it is special compared with other standard things. This folder 
			will be used by r duing the installation. To my understanding,
			this folder will not be CHECKED and subdirectories in this
			folder will be copied to the target location during installation
			and moved up by one level to the top-level directory. This means
			we basically can put things or folders that we want to be deployed
			to the target installation location in here. One thing, however, needs
			attention is that things insides this can not have the names like
			"R" or "DISCRIPTION",etc. 
			Put raw data in here ./inst/extdata!!!!
	./dev -- this folder will not be deployed and it is for developing purpose.
			Code to maintain the package and update dataset (or others ?) is
			kept in here. Also, we will put references, notes, etc in here too.
				notes.txt, more about thoughts of the project
				build_package_to_delete.txt, code to run to install, maintain
				update, build, install the package.
	There also other components and will be added later on.
			