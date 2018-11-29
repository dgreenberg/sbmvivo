#important: try not to edit this file in the matlab editor!
#the matlab editor inserts tabs as spaces which can mess up the makefile

# Download CUDA unbound (cub) and adjust CUBINC variable

CUDAINC := /usr/include
CUDALIB := /usr/lib/x86_64-linux-gnu

CUPFDIR := ../../+sbm/cudasrc
CUBINC  := /data/share/cub
# CUBINC  := /data/src/cub-1.5.2

INCDIRS := -I$(CUDAINC) -I$(CUPFDIR) -I$(CUBINC)
LIBDIRS	:= -L$(CUDALIB)
RPATHS  := -rpath,$(CUDALIB)                        #for use with -Wl
	
# Mex script installed by Matlab, you may have to modify the path
# for Mex to work correctly, you will have to modify mexopts.sh to use the correct versions of gcc/g++
# typically mexopts.sh is found in /home/<USER>/.matlab/<VERSION>/mexopts.sh
# for example, /home/greenberg/.matlab/R2014b/mexopts.sh
# If not present, create via 
#	mex -setup -f /home/<USER>/.matlab/<VERSION>/mexopts.sh
#Probably not necessary if system default is changed as described below
MEX = /usr/local/bin/mex

#this assumes that symbolic links to the correct version of gcc/g++ (should be 4.7 as of Matlab 2014b) have been placed in $(CUDA)/bin
#NOTE: Seems to use system default, not $(CUDA)/bin. System default can be set to gcc-4.7 for the first time by the following commands:
#	sudo apt-get install gcc-4.7 g++-4.7
#	sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6 10
#	sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.7 20
#	sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-6 10  
#	sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.7 20
#CAUTION: It is verz important to reverse this after compilation by
#	sudo update-alternatives --config gcc
#	sudo update-alternatives --config g++
# If not reversed, kernel updates maz fail and machine may not boot properly
#
# NOTE: Up to at least 2017b, MATLAB still only supports gcc4, which is not included in Ubuntu 17.10+. 
# Version 4.8 can be obtained by adding the toolchain/test repository. This tested working with MATLAB 2017a.
# sudo add-apt-repository ppa:ubuntu-toolchain-r/test
NVCCFLAGS := --compiler-bindir=/usr/bin -DNDEBUG -use_fast_math -O=4 -arch=sm_30 --ptxas-options=-v -m 64 -Xcompiler -fPIC -D_FORCE_INLINES
LIBS	  := -lcudart -lcusparse -lcublas

all: 	kernels mex clean

kernels: 
	/usr/bin/nvcc $(CUPFDIR)/cuda5s2b.cu -c -o ./cuda5s2b.cu.o $(INCDIRS) $(NVCCFLAGS)	#use nvcc to compile CUDA code into object files that can be linked from vanilla c++ during mex compilation
	ar -r ./libcuda5s2b.a ./cuda5s2b.cu.o							#create library for static linking 

mex:    
	${MEX} -v LDFLAGS="\$$LDFLAGS -Wall -Wl,$(RPATHS)" -L. -lcuda5s2b cuda5s_mex.cpp $(LIBDIRS) $(LIBS) $(INCDIRS)
	${MEX} -v LDFLAGS="\$$LDFLAGS -Wall -Wl,$(RPATHS)" -L. -lcuda5s2b cuda5s_mex_configure.cpp $(LIBDIRS) $(LIBS) $(INCDIRS)
	${MEX} -v LDFLAGS="\$$LDFLAGS -Wall -Wl,$(RPATHS)" -L. -lcuda5s2b cuda5s_mex_init.cpp $(LIBDIRS) $(LIBS) $(INCDIRS)
	${MEX} -v LDFLAGS="\$$LDFLAGS -Wall -Wl,$(RPATHS)" -L. -lcuda5s2b cuda5s_mex_clear.cpp $(LIBDIRS) $(LIBS) $(INCDIRS)
	${MEX} -v LDFLAGS="\$$LDFLAGS -Wall -Wl,$(RPATHS)" -L. -lcuda5s2b cuda5s_listdevices.cpp $(LIBDIRS) $(LIBS) $(INCDIRS)	
	${MEX} -v LDFLAGS="\$$LDFLAGS -Wall -Wl,$(RPATHS)" -L. -lcuda5s2b cuda5s_devicereset.cpp $(LIBDIRS) $(LIBS) $(INCDIRS)
    
    
clean:
	rm ./libcuda5s2b.a
	rm ./cuda5s2b.cu.o
