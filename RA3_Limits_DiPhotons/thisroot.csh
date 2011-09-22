
unsetenv ROOFITSYS
setenv ROOTSYS /uscmst1b_scratch/lpc1/lpctau/dwjang/work/products/root_v5.30.00
setenv PATH "$ROOTSYS/bin:$PATH"

if ((${?LD_LIBRARY_PATH})) then
       setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}
else
       setenv LD_LIBRARY_PATH ${ROOTSYS}/lib
endif

