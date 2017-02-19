#!/opt/local/bin/python2.7

import pygdml as _pygdml
import numpy as _np



def realistic_collim(vis = "volumes"):

    if (vis=="volumes"):


        """
        PYGDML EXAMPLE - Simple collimator

        All sizes are in implicitly in mm and all angles are in rad. No other units can currently be specified.
        Coordinate system is right handed (rh) with Z along the beam. Positive roation angles produce counter-clockwise (ccw) rotation 

        Solid constructors follow Geant4 primitive solid constructor signatures (not to be mistaken with GDML constrctor signatures).
        For exact solid constructor signatures check constructor docstrings. 

        Example:
            pygdml.Box(name, halflenghtX, halflenghtY, halflenghtZ) #signature of box primitive constructor
            box = pygdml.solid.Box("abox", 7,5,3)                   #box solid object

        A PYGDML volume in the geometry model is represented by a solid with assigned material and placement (similar to a Geant4 physical volume).
        The volume constructor has the following signature:

        pygdml.Volume(self, rot, tlate, currentVolume, name, motherVolume, copyNr, surfCheck, material, scale=[1, 1, 1])

        Where:
           rot           = (list of floats) Rotation about the xyz axes, ex. rot = [0,0,pi/2] for a 90 degree rotation about the Z axis
           tlate         = (list of floats) Postion of the volume in the xyz plane, ex. tlate=[10,0,3] to place a volume at X=10, Y=0, Z=3
           currentVolume = (solid object)   Solid being placed, ex. currentVolume = BOX where BOX is an object of type pygdml.solid.Box
           motherVolume  = (volume object)  The parent volume of the current volume, exmotherVolume = BOX_VOL where BOX_VOL is an object of type pygdml.Volume. 
                                             NB: For top level container (world_volume) montherVoume = None   
           copyNr        = (int)            Integer that keeps track of numbers of identical volumes placed. Set copyNr = 1 as it is present just of parity with Geant4 and is not currently used.
           surfCheck     = (bool)           Surface check flag. Set surfCheck = False as it is present just of parity with Geant4 and is not currently used.
           material      = (string)         Material assigned to volume. Only stored in PYGDML for the purposes of writing out to the output GDML file and nonexistent materials prodce no warnings. ex. material="cheese"
           scale         = (list of floats) Scaling factor for xyz dimensions of the volume. Default is [1,1,1] so no need to set unless needed. 

        The geometry heirarchy is set at construction time and always begins with a "world volume" which has no parent volume. Meshing of volumes for rendering and is invoked through the pycsgmesh() functions

        """
        
        ###SOLID DEFINITIONS AND VOLUME HIEARCHY BUILDING###                
        # Jaw solid is a trapezoid parametrized in terms of half-lengths of faces
        
        l_x = 24     #length of middle box along X
        l_y = 24     #length of middle box along Y
        l_z = 100    #length of middle box along Z

        pDz = 69/2
        pDx1  = 12
        pDx2  = pDx1
        pDx3  = pDx1
        pDx4  = pDx1
        pDy1  = 10/2
        pDy2  = 24/2
        pTheta = _np.arctan2(pDy2,(2*pDz)) 
        pPhi  = 0 
        pAlp1 = 90
        pAlp2 = 90
       
        h     = 8   #Gap between upper jaw and beam pipe center along Y
        gap   = 24  #Separation of collimator jaw center along Y
        bp_r  = 5   #beampipe radius


        """
        Front view of trapezoid jaw piece:

               l_xneg
            +++++++++++++
            +           + l_yneg
            +           + 
            +           +
            +           +       
            +++++++++++++

       Side view of trapezoid jaw piece:

                      
                       ++
                    +   +
                 +      + 
              +         +
          +             +       
       +                +
l_yneg +                + l_ypos         
       ++++++++++++++++++     
          l_z_trd


       """

        #World volume is created first to initalise the volume hierarcy
        world_solid    = _pygdml.solid.Box("world_solid",500,500,500)
        world_volume   = _pygdml.volume.Volume(None, [0,0,0],world_solid, "world_volume", None, 0, False, "G4_Galactic")

        #solid.Tubs(const string& pName, double pRMin, double pRMax, double pHalfLength_z, double pStartPhi, double pSegmentPhi) 
        bp_solid       = _pygdml.solid.Tubs("bp_solid",0, bp_r,(bp_r*1.1),0,2*_np.pi)                              #Makes a cylindrical cutout for the beampipe in the container volume so the colimator can be placed
        world_sub      = _pygdml.solid.Subtraction("world_sub",world_solid,bp_solid,[[0,0,0],[0,0,0]])                   #around the beampipe without overlaps in Geant4
        world_volume   = _pygdml.volume.Volume(None, [0,0,0],world_sub, "world_volume", None, 0, False, "G4_Galactic")   #World volume has parent volume set to None
        
        #Solid definition
        #solid.Box(const string& pName, double pHalfLength_X, double pHalfLength_Y, double pHalfLength_Z) 
        jaw1_solid_box      = _pygdml.solid.Box("jaw_solid_box", l_x/2, l_y/2, l_z/2)
        #solid.Trd(const string& pName, double pHalfLength_X1 (at -dz), double pHalfLength_X2 (at +dz), double pHalfLength_Y1 (at -dz), double pHalfLength_Y2 (at +dz), double pHalfLength_Z (along z)) 
        #solid.Trap(const string& pName, double pHalfLength_Z, double pTheta (polar angle in YZ plane), double pPhi (angle in XZ plane), double pHalfLength_Y1 (at -dz), double pHalfLength_X1 (at -dy1 and -dz), double pHalfLength_X2 (at +dy1 and -dz), double pAlpha1 (Angle between y axis and centre of lower z endcap face), double pHalfLength_Y2 (at +dz), double pHalfLength_X3 (at -dy2 abd +dz), double pHalfLength_X4 (at +dy2 and +dz), double pAlpha2 (Angle between y axis and centre of upper z endcap face) ) 
        jaw1_solid_trd1     = _pygdml.solid.Trap("jaw_solid_trd1", pDz, pTheta, pPhi, pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3, pDx4, pAlp2)
        
        #Union(Name,object1,object2,[[rot2],[transl2]])
        jaw1_solid_trd = _pygdml.solid.Union("jaw_solid_trd",jaw1_solid_trd1, jaw1_solid_trd1,[[0,_np.pi,0],[0,0,pDz*2+l_z]])
        jaw1_solid     = _pygdml.solid.Union("bot_jaw_solid",jaw1_solid_trd, jaw1_solid_box,[[0,0,0],[0,0,pDz+l_z/2]])
        
        #jaw2_solid     = jaw1_solid

        jaw_top_rot    = [0,0,0]
        #jaw_bot_rot    = [0,0,-_np.pi]

        jaw_top_pos    = [0,h+l_y/2,0]
        #jaw_bot_pos    = [0,(h+l_y/2)-gap,0]

        
        jaw_top_volume = _pygdml.volume.Volume(jaw_top_rot, jaw_top_pos, jaw1_solid, "jaw_top", world_volume, 0, False, "G4_Fe") #Volume defintion = solid+placement+hierarcy
        #jaw_bot_volume = _pygdml.volume.Volume(jaw_bot_rot, jaw_bot_pos, jaw2_solid, "jaw_bot", world_volume, 0, False, "G4_Fe") 

        ext = world_volume.setClip()         #Trims the world container to the extent of the scene
        print "\nScene extent: ", ext,"\n"

        m = world_volume.pycsgmesh()         #Meshing method called on world volume invokes recursive meshing for the whole scene to enable visualisation rendering.
                                             #Returns list of all meshes for daugher volumes.
        
        v = _pygdml.VtkViewer()              #Initialise viewer
        v.addSource(m)                       #Adds a list of meshes to display
        v.view()                             #Diplay scene on screen

        g = _pygdml.Gdml()                   #Initialise GDML writer
        g.add(world_volume)                  #Add volume to the writer
        g.write('./realistic_collim.gdml')   #Write GDML file
        
    else:
        print "Unknown visualisation option: ", vis

        
realistic_collim(vis = "volumes")
