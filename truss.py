import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy
import warnings

class Truss:
  """
  This class represents a truss, with attributes joints and beams. The joints 
  are connected by beams, and each beam has a tension force that acts on 
  joints. There may also be other forces like external force and reaction 
  force acting on the joints. If the truss is in equilibrium, the forces on 
  x and y directions of the joints sum up to zero.
  """
  def __init__(self, jointfile, beamfile):
    """
    This methods initiate a Truss object by setting the attribute joints as
    a list of joint data, in which each element represents on joint, and 
    attribute beams as list of beam data, in which each element represents a 
    beam. In each joint, the element contains index, x position, y position,
    external force and whether it is rigidly supported. In each beam, the 
    element contains index, and the two joints on either side of the beam.
    """
    joints=[]
    beams=[]  
    #read in files  
    try:
      with open(jointfile,'r') as f:
        next(f)
        for line in f:
          joints.append(line.strip().split())
    except RuntimeError:
      print("Error when reading in joints data")
    try:
      with open(beamfile,'r') as f:
        next(f)
        for line in f:
          beams.append(line.strip().split())
    except RuntimeError:
      print("Error when reading in beams data")
    #convert data type and initialize attributes
    self.joints=np.asarray(joints).astype(float)
    self.joints.astype(float)
    self.beams=np.asarray(beams).astype(float).astype(int)
  
  def calculate_force(self):
    """
    This method calculates the tension force on each beam if applicable
    and returns a list of such forces in order.
    """
    #create csr_matrix components
    row=[]
    col=[]
    data=[]
    #add beam coefficients
    for b in range(self.beams.shape[0]):
      #find two joints of the beam
      j1=self.joints[self.beams[b][1]-1]
      j2=self.joints[self.beams[b][2]-1]
      j1num=int(j1[0])
      j2num=int(j2[0])
      #calculate b coefficients on joint 1
      j1x=j1[1]-j2[1]
      j1y=j1[2]-j2[2]
      bj1x=j1x/np.sqrt(j1x**2+j1y**2)
      bj1y=j1y/np.sqrt(j1x**2+j1y**2)
      #add data to sparse matrix
      row.append(j1num*2-2)
      col.append(b)
      data.append(bj1x)
      row.append(j1num*2-1)
      col.append(b)
      data.append(bj1y)
      #calculate b coefficients on joint 2
      j2x=j2[1]-j1[1]
      j2y=j2[2]-j1[2]
      bj2x=j2x/np.sqrt(j2x**2+j2y**2)
      bj2y=j2y/np.sqrt(j2x**2+j2y**2)
      #add data to sparse matrix
      row.append(j2num*2-2)
      col.append(b)
      data.append(bj2x)
      row.append(j2num*2-1)
      col.append(b)
      data.append(bj2y)
    #vector on right hand side
    sol=[0]*(self.joints.shape[0]*2)
    #add reaction force coefficient and external force
    for j in range(self.joints.shape[0]):
      #add Fx
      if self.joints[j][3]!=0:
        jnum=int(self.joints[j][0])
        sol[jnum*2-2]=-self.joints[j][3]
      #add Fy
      if self.joints[j][4]!=0:
        jnum=int(self.joints[j][0])
        sol[jnum*2-1]=-self.joints[j][4]
      #add reaction force
      if self.joints[j][5]==1:
        jnum=int(self.joints[j][0])
        row.append(jnum*2-2)
        col.append(max(col)+1)
        data.append(1)
        row.append(jnum*2-1)
        col.append(max(col)+1)
        data.append(1)
    beamforce=[]
    #raise warning about overdetermined/underdetermined
    if max(row)>max(col) or max(row)<max(col):
      raise RuntimeError("Truss geometry not suitable for static equilbrium analysis")
    else:
      #create sparse matrix
      matrix=scipy.sparse.csr_matrix((data, (row, col)), shape=(max(row)+1, max(col)+1))
      sol=np.asarray(sol)
      # Catch warnings as exceptions
      warnings.filterwarnings('error')
      try:
          force=scipy.sparse.linalg.spsolve(matrix,sol)
          beamforce=force[:self.beams.shape[0]]  
      except:
          print("Cannot solve the linear system, unstable truss?")
    return beamforce
    
  def PlotGeometry(self, plot_output):
    """
    This method plots the truss based on coordinates of the joints
    and saves the plot to the given output file name.
    """
    for b in range(self.beams.shape[0]):
      #find two joints of the beam
      j1=self.joints[self.beams[b][1]-1]
      j2=self.joints[self.beams[b][2]-1]
      j1x=j1[1]
      j1y=j1[2]
      j2x=j2[1]
      j2y=j2[2]
      x=[j1x,j2x]
      y=[j1y,j2y]
      plt.plot(x,y,color="blue")
      plt.axis('equal')
      plt.margins(0.1)
    plt.savefig(plot_output)
    
  def __repr__(self):
    """
    This method prints the tension force on each beam.
    """
    beamforce=self.calculate_force()
    ret=""
    if len(beamforce)!=0:
      ret=" Beam       Force\n-----------------\n"
      for i in range(len(beamforce)):
        f="{:.3f}".format(beamforce[i])
        if f[0]!='-':
          f=" "+f
        ret+=("    {}{:>12}\n".format(i+1,f))
    return ret
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  