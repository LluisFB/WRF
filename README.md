### WRF-ARW Modeling System  ###

We request that all new users of WRF please register. This allows us to better determine how to support and develop the model. Please register using this form:[https://www2.mmm.ucar.edu/wrf/users/download/wrf-regist.php](https://www2.mmm.ucar.edu/wrf/users/download/wrf-regist.php).

For an overview of the WRF modeling system, along with information regarding downloads, user support, documentation, publications, and additional resources, please see the WRF Model Users' Web Site: [https://www2.mmm.ucar.edu/wrf/users/](https://www2.mmm.ucar.edu/wrf/users/).
 
Information regarding WRF Model citations (including a DOI) can be found here: [https://www2.mmm.ucar.edu/wrf/users/citing_wrf.html](https://www2.mmm.ucar.edu/wrf/users/citing_wrf.html).

The WRF Model is open-source code in the public domain, and its use is unrestricted. The name "WRF", however, is a registered trademark of the University Corporation for Atmospheric Research. The WRF public domain notice and related information may be found here: [https://www2.mmm.ucar.edu/wrf/users/public.html](https://www2.mmm.ucar.edu/wrf/users/public.html).


# Lluís Fita Borrell fork #
This is the WRF repository of Lluis Fita Borrell, forked from WRF v7.1. CIMA-IFAECI, UBA-CONICET-CNRS-IRD, C. A. Buenos Aries Argentina

## Installation ##

Getting the code from the fork in a INSTDIR
```
$ cd $INSTDIR
$ git clone git@github.com:LluisFB/WRF.git
```

Getting submodules
```
$ git submodule update --init --recursive
```

Listing available branches (all official from WRF v7.1) and the ones created by L Fita
```
$ git branch -a
```

To switch to any WRF branch
```
$ git checkout [branch_name]
```
