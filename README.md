## GISAssignment
GIS assignment, visualising the spatial distribution of the UHI effect in Barcelona and Los Angeles using Landsat 8 OLI and TIRS data

Download the whole "Barcelona-GI-Study" folder to run the RStudio script that created all the plots seen in the report. 
You can also access a Rnotebook here:
http://rpubs.com/HReynier_/GISAssignment

# IMPORTANT!!!
As the landsat 8 data is very large, I have been unable to upload it to this github, in which case you must follow the following instructions.

# Barcelona Data Download
Go to https://earthexplorer.usgs.gov/
1. Search for Barcelona in the address/place box and select "Barcelona, Spain".
2. Select the date range between 9/5/2019 and 11/5/2019 - as it is a US website, check dates are correct
for 9th May 2019 - 11th May 2019.
3. Click dataset, select Landsat, Landsat Collection 1 Level-1.
4.Click results and select the one image that is returned. Download "Level-1 GeoTIFF Data Product (869.0 MB)" 
May take a while..
5. Make sure to save the download inside the "Barcelona/May_10_19" folder, and unzip the tar.gz to get 12 tif files.

# Los Angeles Data Download
Go to https://earthexplorer.usgs.gov/
1. Search for Los Angeles in the address/place box and select "Los Angeles, CA, USA".
2. Select the date range between 29/5/2019 and 31/5/2019 - as it is a US website, check dates are correct
for 29th May 2019 - 31st May 2019.
3. Click dataset, select Landsat, Landsat Collection 1 Level-1.
4.Click results and select the one image that is returned. Download "Level-1 GeoTIFF Data Product (869.0 MB)" 
May take a while..
5. Make sure to save the download inside the "LA/May_30_19" folder, and unzip the tar.gz to get 12 tif files.

After this the R Script "Assignment.R" and the RNotebook should run correctly.
