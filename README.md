MECA: MaxEnt Climate Aggregator
===============================

Takes gridded climate models and produces from them the data files necessary to
run the MaxEnt program.

Installation
------------

You will need the following:

 * Python (https://www.python.org/downloads/)
 * h5py   (https://pypi.python.org/pypi/h5py)
 * numpy  (http://www.scipy.org/scipylib/download.html)

On **Linux** these can usually be installed with your package manager using, e.g.,

    sudo apt-get install python-numpy

or

    sudo pip install numpy

And everything will be taken care of.

Binaries for **Windows** and **Mac** should be available through the links
above. It may take some fiddling to get things to work.



Usage
-----

Check [here, at gdo-dcp.ucllnl.org](http://gdo-dcp.ucllnl.org/downscaled_cmip_projections/dcpInterface.html#Projections:%20Complete%20Archives), to get gridded climate files of interest. Also available in [this](ftp://gdo-dcp.ucllnl.org/pub/dcp/archive/cmip5/) FTP repository.

Place the files of interest within a directory. The program will search
recursively through the directory and all of its subdirectories to find climate
files. Because of this, you can keep the files organized in their own folders or
just copy the FTP directory structure from above.

Run the program from the command line using the following

    meca.py <basedir> <rcp> <startyear> <endyear> <output_prefix>

The arguments are:

    <basedir>:       Where you put the gridded climate data
    <rcp>:           This is one of historical/rcp26/rcp45/rcp60/rcp85
    <startyear>:     First year of period of interest
    <endyear>:       Last year of period of interest
    <output_prefix>: Output files will begin with this name/directory