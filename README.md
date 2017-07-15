# reveals-na

This repository contains the code to reconstruct vegetation composition from fossil pollen data for North America using the REVEALS model described in Sugita (2007) and implement in DISQOVER. 

This is a work in progress. The workflow includes:
  - Pulling pollen records for North America from Neotoma
  - Compiling 
  - Aggregating pollen to appropriate time bins
  - Run RevealsInR on all sites (likely want to restrict to only large lakes)
  - Average cells withing each one degree grid cell
  - Map and summarize results
