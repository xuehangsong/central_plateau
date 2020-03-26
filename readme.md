

# Task 2 #

## Object   

**Email from MLR (Thu 3/12/2020 12:03 PM)**
The files were copied to constance located in /pic/projects/dvz/xhs_simus/Hanford_ACMs  

*Will you please take the MODFLOW parameter fields located here \\pnl\projects\HanfordACMs\InteraCentralPlateauCP-47631-Rev4\calibratedmodel\historic_flow
and write out a file containing the following variables for their model grid centroids:
Easting[m], Northing[m], Elevation [m, NAVD88], Ksxx[m/d], Ksyy[m/d], Kszz[m/d], Porosity
You can use charge code NE7119 for this part

We also need this for the Plateau-to-River model located here \\pnl\projects\HanfordACMs\InteraPlateauToRiverCP-57037-Rev2\calibfinal. 
You can use charge code NE8363 for this part.

This will support a nearest neighbor type property assignment. We should also think about and discuss what it would take to consider the two meshes and their volumes of overlap when estimating properties/parameters for one from the other. 
*

## highlights  

**1. model discretization in horizontal direction**  
The horizontal direction discretization is read from *.spc

"*.dis" is the discretization file for MODFLOW 2000, while I didn't find x,y origin (XOFF, YOFF) in this file  
"* .spc" is the  grid specification file(?). It has all the discretization information.

**2. model layers**  
The names of model layer files (top, bot files) was consistent with the file list in  "*.dis"  

**3. model ss, sy, hz,hy**  
This part has to be hard coded. 
the list of ss,sy,hz,hy is hard wired based the information in "*lpf" file

**4. I assume porosity equals to sy**  

