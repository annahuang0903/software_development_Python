import glob
import math
import os
"""
This class is designed to represent airfoils with the Cartesian coordinates of 
panels, and calculate the lift coefficients and stagnation points given the 
pressure coefficient data with different angles of attack. For each Airfoil 
object, the attributes include the directory of all relevant data, xy coordinates
of all points defining the panels, the chord length (which is the distance 
from the leading edge to the trailing edge), the list of alphas given in different 
pressure coefficient files, and the list of pressure coefficient file names. 
"""

class Airfoil:
  def __init__(self,inputdir):
    """
    This function initializes an Airfoil object with a given data directory 
    and call the load_data() function to initialize all its attributes.
    """
    #if inputdir does not exist in current directory
    if len(glob.glob(inputdir))==0:
      raise RuntimeError("Not a valid directory")
    self.load_data(inputdir)
    
  def load_data(self,inputdir):
    """
    This function reads in the xy coordinates of panel points, calculates the 
    chord length, reads in the list of file names containing the pressure 
    coefficients, and the list of alphas corresponding to those files.
    """
    #xy read in
    self.inputdir=inputdir
    cur_path = os.path.dirname(__file__)
    xyfile_path=cur_path+"/"+self.inputdir+"xy.dat"
    xyfile=[]
    #handling error when reading in file
    try:
      with open(xyfile_path,'r') as f:
        next(f)
        for line in f:
          xyfile.append(line.strip().split())
    except IOError:
      print("Error when reading in XY data")
    xy=[]
    for i in range(len(xyfile)):
      xy.append([float(j) for j in xyfile[i]])
    self.x=[]
    self.y=[]
    for item in xy:
      self.x.append(item[0])
      self.y.append(item[1])
    #chord length
    self.chord=max(self.x)-min(self.x)
    #read in alphas
    filelist=os.listdir(cur_path+"/"+inputdir)
    filelist.sort()
    alphafile=filelist[:-1]
    self.alphas=[]
    for i in alphafile:
      self.alphas.append(float(i[5:][:-4]))
    #sort alpha files based on alpha
    self.alphas.sort()
    self.alphas_sign=[]
    for i in self.alphas:
      if i<0:
        self.alphas_sign.append(str(i))
      elif i==0:
        self.alphas_sign.append(str(i))
      else:
        alpha_name='+'+str(i)
        self.alphas_sign.append(alpha_name)
    #get all the file names containing alpha
    self.alphafile=[]
    for i in range(len(self.alphas_sign)):
      file=alphafile[0][:5]+self.alphas_sign[i]+alphafile[0][-4:]
      self.alphafile.append(file) 
    
  def read_cp(self,alphafile):
    """
    This function reads in the particular pressure coefficient file and returns
    the list of pressure coefficients in that file.
    """
    cp=[]
    cur_path = os.path.dirname(__file__)
    alphafile_path=cur_path+"/"+self.inputdir+alphafile
    #handling error when reading in file
    try:
      with open(alphafile_path,'r') as f:
        next(f)
        for line in f:
          cp.append(line.strip())
    except IOError:
      print("Error when reading in alpha data")
    cp=[float(i) for i in cp]
    return cp
  
  def calc_cl(self,cp,alpha):
    """
    This function calculates the lift coefficients of the airfoil under a given 
    angle of attacks. It returns the lift coefficient.
    """
    cx=0
    cy=0
    for i in range(len(cp)):
      dx=-(self.y[i+1]-self.y[i])*cp[i]/self.chord
      cx+=dx
      dy=(self.x[i+1]-self.x[i])*cp[i]/self.chord
      cy+=dy
    cl=cy*math.cos(alpha/180*math.pi)-cx*math.sin(alpha/180*math.pi)
    return cl
  
  def stag_pt(self,cp):
    """
    This function finds the stagnation point of the air foil under a specific 
    angle of attack. The list of pressure coefficients is passed in and the
    coordinates of the stagnation point and the pressure coefficient closest 
    to 1 are returned.
    """
    max_cp=max(cp)
    cp_max=cp.index(max_cp)
    x_stag=(self.x[cp_max]+self.x[cp_max+1])/2  
    y_stag=(self.y[cp_max]+self.y[cp_max+1])/2 
    return x_stag,y_stag,max_cp
  
  def __repr__(self):
    """
    This function prints all the calculated lift coefficients under all angles 
    of attack given in the file directory, along with the stagnation points
    and maximum value of pressure coefficients.
    """
    string="Test case: NACA {}\n\n".format(self.inputdir[4:8])
    string+="  alpha     cl          stagnation pt\n"
    string+="  -----  -------  --------------------------\n"
    for i in range(len(self.alphas)):
      cp=self.read_cp(self.alphafile[i])
      cl=self.calc_cl(cp,self.alphas[i])
      if cl>=0:
        cl=" {:.4f}".format(cl)
      else:
        cl="{:.4f}".format(cl)
      if self.alphas[i]>=0:
        alpha=" "+str(self.alphas[i])
      else:
        alpha=str(self.alphas[i])
      x_stag,y_stag,max_cp=self.stag_pt(cp)
      if x_stag>=0:
        x_stag=" {:.4f}".format(x_stag)
      else:
        x_stag="{:.4f}".format(x_stag)
      if y_stag>=0:
        y_stag=" {:.4f}".format(y_stag)
      else:
        y_stag="{:.4f}".format(y_stag)
      string+="  {}   {}  ({}, {})  {:.4f}\n".format(alpha,cl,x_stag,y_stag,max_cp)
    return string
 
