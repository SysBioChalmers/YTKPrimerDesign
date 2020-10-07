# YTKPrimerDesign
This repo contains two shiny apps for primer design any an outdated python app.

To run the shiny apps locally on your computer, do the following:
1. Make sure you have R and RStudio installed
2. Install the shiny package:
install.packages("shiny")
3. Download this repo somewhere on your disk.
4. Open the app.R file for your app (there are two apps in different folders)
5. Press the button "Run App" in the top right corner of the editor window.




The python program is a bit outdated and is not recommended, it was intended to do the same as one of the apps.
It requires the python package primer3-py.

To run the python program in PyCharm, do the following:

1. Clone the repo to disk, install python and PyCharm.
2. Open PyCharm, select file->open in the menu and choose the directory PythonGui
3. Configure an interpreter. This is done by selecting file->settings. Then go to
   project: PythonGui, and under that Project interpreter. In the combo box, you
   can select an interpreter. Click show all there. A new dialog opens. Click + and 
   add an interpreter it suggests, typically Python 3.6 or something similar. Then
   make sure that it is selected in the combobox when you get back to the previous
   window.
4. On the same settings page, you can install packages. You need to install primer3-py
   by clicking the + sign next to the list. The primer3 package should be available 
   in the list, so just install it.
5. It should now be possible to run the program by selecting Run->Run MoClo
