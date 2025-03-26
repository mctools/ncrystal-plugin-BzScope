# BzScope Plugin for NCrystal

This repository contains code for the NCrystal plugin named `BzScope`. This plugin is able to load the coherent phonon scattering scattering functions that calculated by the tool BzScope.

The main development tree at https://code.ihep.ac.cn/cinema-developers/ncplugin-bzscope.git


## Installation

### Prerequisites
This plugin requires the following NCrystal packages:
```bash
pip install ncrystal ncrystal-python ncrystal-pluginmanager
```


### Install the Plugin
```bash
pip install git+https://code.ihep.ac.cn/cinema-developers/ncplugin-bzscope
```

### Usage Example
To test the plugin with a sample BeO material file:
```bash
nctool plugins::BzScope/bzscope_beo_c1_300K.ncmat
```
This will display the results from the example BeO data file.



It is currently under development and not ready for general usage.
