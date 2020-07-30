# pines_analysis_toolkit

[PINES website](https://pines.bu.edu/)

This is a package for working with PINES data.

INSTALLATION FOR GENERAL USE. 
  1. Uninstall the current version of pines_analysis_toolkit on your computer, if you have one installed: 
        pip (or pip3) uninstall pines_analysis_toolkit
     Type 'y' when prompted to complete the uninstall. 
     
  2. Install the new version of the pipeline from github, using: 
        pip (or pip3) install --upgrade git+git://github.com/patricktamburo/pines_analysis_toolkit.git
  
  3. Launch python, and see if you can import pines_analysis_toolkit: 
        python (or python3)
        
        import pines_analysis_toolkit as pat
        
     If the import works without throwing any errors, the pipeline has installed successfully!
   
   4. Check your ~/Documents/ folder for a new directory called PINES_analysis_toolkit/. All PINES data products will be placed in this directory. 
  
INSTRUCTIONS FOR CONTRIBUTING TO CODEBASE.
  1. Having completed the above, next create a directory where you will work on the pipeline. Mine is in ~/Documents/python_programs/pines_analysis_toolkit_package/. When writing new features for the toolkit, you will do so in this directory.
  
  2. Navigate a terminal to that directory, and type the following: 
  
        git init
        
        git remote add origin https://github.com/patricktamburo/pines_analysis_toolkit.git
        
        git pull origin master
        
      Check to make sure files have appeared in the directory.
    
  3. Start python from your equivalent of ~/Documents/python_programs/pines_analysis_toolkit_package/, and type import pines_analysis_toolkit as pat. If it imports without throwing any errors, you're all set!
  
  4. Create a new branch to work on whatever feature you're trying to implement. From terminal: 
    
      git checkout -b "new-branch"
      
  5. Write/test your new feature.
  
  6. When the feature is performing to your liking, make a new commit of your changes: 
      git add . 
      
      git commit -m "added new feature, fixed bug, etc."
      
  7. Push your edited branch to github: 
      git push remote new-branch
      
  8. Once your changes are merged into the master branch, you can follow the INSTALLATION FOR GENERAL USE instructions to reinstall the package with the updated feature. In your ~/Documents/python_programs/pines_analysis_toolkit_package/ directory, you can git pull origin master to make sure everything is up to date, then repeat steps 3-8 until we discover some Earth-analogs around brown dwarfs. 
  

MODULES (all still under development)
 
    data - programs to grab raw data, science data, calibrations, logs, reducing data, etc. 
    
    photometry - programs for performing photometry on reduced data
  
    analysis - programs for making lightcurves, analyzing variability, etc.
    
    utils - programs for various quality-of-life things
