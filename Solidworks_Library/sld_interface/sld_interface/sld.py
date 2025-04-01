import subprocess as sb
import win32com.client
import pythoncom
import numpy as np
import os
import array
import math

def startSW():
    ## Starts Solidworks
    SW_PROCESS_NAME = r'C:/Program Files/SOLIDWORKS Corp/SOLIDWORKS/SLDWORKS.exe'
    sb.Popen(SW_PROCESS_NAME)

def shutSW():
    ## Kills Solidworks
    sb.call('Taskkill /IM SLDWORKS.exe /F')

def connectToSW():
    ## With Solidworks window open, connects to application      
    sw = win32com.client.Dispatch("SLDWORKS.Application")
    return sw

def sw_to_step_file(sw_file_path, step_file_path):
    global sw
    sw.Visible = False
    error = win32com.client.VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    warning = win32com.client.VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    doc = sw.OpenDoc6(sw_file_path, 1, 0, "", error, warning)
    activeDoc = sw.ActiveDoc
    exportData = sw.GetExportFileData(1)
    success = activeDoc.SaveAs3(step_file_path, 0, 2)
    sw.CloseDoc(activeDoc.GetTitle)

def sw_to_ps_file(sw_file_path, ps_file_path):
    global sw
    sw.Visible = False
    error = win32com.client.VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    warning = win32com.client.VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    doc = sw.OpenDoc6(sw_file_path, 1, 0, "", error, warning)
    activeDoc = sw.ActiveDoc
    exportData = sw.GetExportFileData(2)
    success = activeDoc.SaveAs3(ps_file_path, 0, 2)
    sw.CloseDoc(activeDoc.GetTitle)    

def ps_to_stl(file_name, stl_file_path):
    global sw, swPartDoc
    error = win32com.client.VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    warning = win32com.client.VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    doc = sw.LoadFile2(file_name, "")    
    activeDoc = sw.ActiveDoc
    body_array = activeDoc.GetBodies2(0,True)
    body = body_array[0]
    mesh_body = body.ConvertToMeshBody(False, False, False, 0.005, True, 0.01,0.5,True, 0.0002,-1)
    success = activeDoc.SaveAs3(stl_file_path, 0, 2)
    sw.CloseDoc(activeDoc.GetTitle)

def ps_to_stl_single_helix(file_name, stl_file_path):
    global sw, swPartDoc
    error = win32com.client.VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    warning = win32com.client.VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    doc = sw.LoadFile2(file_name, "")    
    activeDoc = sw.ActiveDoc
    body_array = activeDoc.GetBodies2(0,True)
    body = body_array[0]
    mesh_body = body.ConvertToMeshBody(False, False, False, 0.005, True, 0.0000073,0.472,True, 0.00012,-1)
    success = activeDoc.SaveAs3(stl_file_path, 0, 2)
    sw.CloseDoc(activeDoc.GetTitle)


def ps_to_stl_ir2(file_name, stl_file_path):
    global sw, swPartDoc
    error = win32com.client.VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    warning = win32com.client.VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    doc = sw.LoadFile2(file_name, "")    
    activeDoc = sw.ActiveDoc
    body_array = activeDoc.GetBodies2(0,True)
    body = body_array[0]
    mesh_body = body.ConvertToMeshBody(False, False, False, 0.005, True, 0.01,0.5,True, 0.00011,-1)
    success = activeDoc.SaveAs3(stl_file_path, 0, 2)
    sw.CloseDoc(activeDoc.GetTitle)

def openFile(sw, Path):
    ## With connection established (sw), opens part, assembly, or drawing file            
    f = sw.getopendocspec(Path)
    model = sw.opendoc7(f)
    return model

def updatePrt(model):
    ## Rebuilds the active part, assembly, or drawing (model)
    model.EditRebuild3

def zoom_fit():
    global swPartDoc
    swPartDoc.ViewZoomtofit2()

def trimetric_zoom_fit():
    global swPartDoc
    swPartDoc.ShowNamedView2("*Trimetric", -1)
    swPartDoc.ViewZoomtofit2()
  
def takeeN(No):
    return No

def mean(arr):
    return np.mean(arr)

def getSW():
    return sw

def exit_sw():
    global sw
    sw.ExitApp()
    return "App Exited"

def save_step_file(save_location):
    global swPartDoc, sw
    #sw.Visible = True
    activeDoc = sw.ActiveDoc
    exportData = sw.GetExportFileData(1)
    success = activeDoc.SaveAs3(save_location, 0, 2)

def save_sw_file(save_location):
    global swPartDoc, sw
    swPartDoc.SaveAs3(save_location,0,0)
    return "File Saved"

def move(previous_row, axis, rel_y_shift,rel_theta_shift):
    global swPartDoc, sw
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByID2(previous_row, "SOLIDBODY", 0, 0, 0, False, 1, pythoncom.Nothing, 0)
    translated_row = swPartDoc.FeatureManager.InsertMoveCopyBody2(0, rel_y_shift, 0, 0, 0, 0, 0, 0, 0, 0, False, 1)    
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByID2(translated_row.name, "SOLIDBODY", 0, 0, 0, False, 1, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2(axis, "AXIS", 0, 0, 0, True, 2, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2(axis, "AXIS", 0,0,0, True, 4, pythoncom.Nothing, 0)
    new_row = swPartDoc.FeatureManager.InsertMoveCopyBody2(0, 0, 0, 0, 0, 0, -0.0015, 0, rel_theta_shift, 0, False, 1)
    return new_row.name

def copy_move_v2(previous_row, axis, rel_y_shift,rel_theta_shift):# This V2 version is to allow just theta shift without y translation of completely built HS1 helix. 
    global swPartDoc, sw
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByID2(previous_row, "SOLIDBODY", 0, 0, 0, False, 1, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2(axis, "AXIS", 0, 0, 0, True, 2, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2(axis, "AXIS", 0,0,0, True, 4, pythoncom.Nothing, 0)
    new_row = swPartDoc.FeatureManager.InsertMoveCopyBody2(0, 0, 0, 0, 0, 0, -0.0015, 0, rel_theta_shift, 0, True, 1)# -0.0015 is the 1.5 mm centre point about which we are rotating. we need to make this variable is the relative location of the sketch and cylinder centre are changd in any way. 
    return new_row.name

def copy_move(previous_row, axis, rel_y_shift,rel_theta_shift):
    global swPartDoc, sw
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByID2(previous_row, "SOLIDBODY", 0, 0, 0, False, 1, pythoncom.Nothing, 0)
    translated_row = swPartDoc.FeatureManager.InsertMoveCopyBody2(0, rel_y_shift, 0, 0, 0, 0, 0, 0, 0, 0, True, 1)
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByID2(translated_row.name, "SOLIDBODY", 0, 0, 0, False, 1, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2(axis, "AXIS", 0, 0, 0, True, 2, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2(axis, "AXIS", 0,0,0, True, 4, pythoncom.Nothing, 0)
    new_row = swPartDoc.FeatureManager.InsertMoveCopyBody2(0, 0, 0, 0, 0, 0, -0.0015, 0, rel_theta_shift, 0, False, 1)
    return new_row.name

def copy_move_vro(previous_row, axis, rel_z_shift,rel_theta_shift):# This version of copy_move is for resolute HS2 stent main helix
    global swPartDoc, sw
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByID2(previous_row, "SOLIDBODY", 0, 0, 0, False, 1, pythoncom.Nothing, 0)
    translated_row = swPartDoc.FeatureManager.InsertMoveCopyBody2(0, 0, rel_z_shift, 0, 0, 0, 0, 0, 0, 0, True, 1)
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByID2(translated_row.name, "SOLIDBODY", 0, 0, 0, False, 1, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2(axis, "AXIS", 0, 0, 0, True, 2, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2(axis, "AXIS", 0,0,0, True, 4, pythoncom.Nothing, 0)
    new_row = swPartDoc.FeatureManager.InsertMoveCopyBody2(0, 0, 0, 0, 0, 0, 0, rel_theta_shift, 0, 0, False, 1)
    return new_row.name

def reorient_ro(stent_name):# Reorient single helix stent to -y axis
    global swPartDoc, sw
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByID2(stent_name, "SOLIDBODY", 0, 0, 0, False, 1, pythoncom.Nothing, 0)
    rotated_stent = swPartDoc.FeatureManager.InsertMoveCopyBody2(0, 0, 0, 0, 0, 0, 0, 0, 0, -1.5707963267949, False, 0)
    return rotated_stent.name

def row_combine(previous_row, current_row):
    global swPartDoc
    swPartDoc.Extension.SelectByID2(previous_row, "SOLIDBODY", 0, 0, 0, False, 2, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2(current_row, "SOLIDBODY", 0, 0, 0, True, 2, pythoncom.Nothing, 0)
    new_row = swPartDoc.FeatureManager.InsertCombineFeature(15903, pythoncom.Nothing, pythoncom.Nothing)
    return new_row.name

def circpattern(stent,axis,N_rep):
    global swPartDoc
    swPartDoc.Extension.SelectByID2(stent, "SOLIDBODY", 0, 0, 0, False, 256, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2(axis, "AXIS", 0, 0, 0, True, 1, pythoncom.Nothing, 0)
    swCpData = swPartDoc.FeatureManager.CreateDefinition(5)
    swCpData.Direction2 = False
    swCpData.EqualSpacing = True
    swCpData.GeometryPattern = False
    swCpData.ReverseDirection = False
    swCpData.TotalInstances = N_rep
    swCpData.VarySketch = False
    swCpData.spacing = 2*math.pi
    circ_pattern = swPartDoc.FeatureManager.CreateFeature(swCpData)
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByID2(stent, "SOLIDBODY", 0, 0, 0, False, 2, pythoncom.Nothing, 0)
    if N_rep == 2:
        body_name = circ_pattern.name
        swPartDoc.Extension.SelectByID2(body_name, "SOLIDBODY", 0, 0, 0, True, 2, pythoncom.Nothing, 0)
    else:
        for j in range(int(N_rep-1)):
            body_name = circ_pattern.name + "[" + str(j + 1) + "]"
            swPartDoc.Extension.SelectByID2(body_name, "SOLIDBODY", 0, 0, 0, True, 2, pythoncom.Nothing, 0)
    
    row = swPartDoc.FeatureManager.InsertCombineFeature(15903, pythoncom.Nothing, pythoncom.Nothing)
    return row.name

def create_circle(circle_data):
    global swSketchManager
    swSketchManager.CreateCircle(0, 0, 0, circle_data[1][1], 0, 0)
    return str(type(circle_data))

def create_stent_segment_width(segID, width):
    global swPartDoc, swSketchManager
    swPartDoc.Extension.SelectByID2(segID, "SKETCHSEGMENT", 0, 0, 0, False, 0, pythoncom.Nothing, 0)
    swSketchManager.SketchOffset2(width/2, True, True, 2, 1, False)
    return "Created Width"

def combine_subtract(cylinder_name,wrap_name):
    global swPartDoc
    swPartDoc.Extension.SelectByID2(cylinder_name, "SOLIDBODY", 0, 0, 0, False, 1, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2(wrap_name, "SOLIDBODY", 0, 0, 0, True, 2, pythoncom.Nothing, 0)
    stent_one = swPartDoc.FeatureManager.InsertCombineFeature(15902, pythoncom.Nothing, pythoncom.Nothing)
    return stent_one.Name

def create_base_cylinder_axis():
    global swPartDoc
    swPartDoc.Extension.SelectByRay(0, 0, 0, 0, 0, -1, 0.1, 2, False, 1, 0)
    swPartDoc.InsertAxis2(True)

def create_z_axis():
    global swPartDoc
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0, 0, 0, True, 0, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2("Front Plane", "PLANE", 0, 0, 0, True, 0, pythoncom.Nothing, 0)
    status = swPartDoc.InsertAxis2(True)
    return status

def wrap_sketch_around_base_cylinder(SketchName, stent_thickness):
    global swPartDoc, swSketchManager
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByRay(0, 0, 0, 0, 0, -1, 0.1, 2, False, 1, 0)
    swPartDoc.Extension.SelectByID2(SketchName, "SKETCH", 0, 0, 0, True, 4, pythoncom.Nothing, 0)
    wrap = swPartDoc.FeatureManager.InsertWrapFeature2(1, stent_thickness, False, 0, 5)
    return wrap.Name

def wrap_sketch_around_base_cylinder_v2(SketchName, stent_thickness):
    global swPartDoc, swSketchManager
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByRay(0, 0, 0, 0, 0, -1, 0.1, 2, False, 1, 0)
    swPartDoc.Extension.SelectByID2(SketchName, "SKETCH", 0, 0, 0, True, 4, pythoncom.Nothing, 0)
    wrap = swPartDoc.FeatureManager.InsertWrapFeature2(1, stent_thickness, False, 0, 5)
    faces = wrap.GetFaces
    if faces is None:
        return "NOWRAP"
    else:
        return wrap.Name

def wrap_sketch_around_base_cylinder_v3(SketchName, stent_thickness):# Wrapping specifically for bottom strut of HS1 as selection byRay method cannot be further used. 
    global swPartDoc, swSketchManager
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByID2(SketchName, "SKETCH", 0, 0, 0, False, 4, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByRay(0, 1.5E-3, 0, 0, 0, -1, 0.1, 2, True, 1, 0)# Shifted the ray starting point to select the cylinder face instead of previous wrap 
    wrap = swPartDoc.FeatureManager.InsertWrapFeature2(1, stent_thickness, False, 0, 5)
    faces = wrap.GetFaces
    if faces is None:
        return "NOWRAP"
    else:
        return wrap.Name

def delete_sketch(SketchName):
    global swPartDoc, swSketchManager
    swPartDoc.Extension.SelectByID2(SketchName, "SKETCH", 0, 0, 0, True, 4, pythoncom.Nothing, 0)
    swPartDoc.EditDelete()
    return "Sketch Deleted"

def create_base_cylinder(radius, stent_length):
    global swPartDoc, swSketchManager
    topPlane = swPartDoc.FeatureByName("Top Plane")
    topPlane.Select2(False, -1)
    swSketchManager.InsertSketch(True)# Insert Sketch
    swSketchManager.AddToDB = True# These two lines added to allow the circle less than 1 mm radius to be made
    swSketchManager.DisplayWhenAdded = False
    base_circle = swSketchManager.CreateCircleByRadius(0, 0.0015, 0, radius)# Create a circle on the plane at y = 1.5 mm with the radius as received
    swSketchManager.AddToDB = False
    swSketchManager.DisplayWhenAdded = True
    swSketchManager.InsertSketch(True)# exit sketch
    swPartDoc.Extension.SelectByID2(base_circle.GetName, "SKETCHSEGMENT", 0, 0, 0, False, 0, pythoncom.Nothing, 0)# Select circle sketch
    cylinder = swPartDoc.FeatureManager.FeatureExtrusion2(False, False, False, 0, 0, 0.2*stent_length, 1.2*stent_length, False, False, False, False, 0, 0, False, False, False, False, False, True, True, 0, 0, False)
    return cylinder.Name

def create_spline(data_points):
    global swSketchManager
    arg1 = win32com.client.VARIANT(pythoncom.VT_ARRAY | pythoncom.VT_R8, data_points)
    Nothing = pythoncom.Nothing
    argout = win32com.client.VARIANT(pythoncom.VT_BYREF | pythoncom.VT_VARIANT, 0)
    spline = swSketchManager.CreateSpline3(arg1,Nothing,Nothing,True,argout)
    return spline.GetName

def create_spline_3d(data_points, prev_spline_name):
    global swSketchManager, swPartDoc
    i = 0
    swPartDoc.InsertCurveFileBegin()
    for j in range(int(len(data_points)/3)):
        swPartDoc.InsertCurveFilePoint(data_points[i], data_points[i + 1], data_points[i + 2])
        i = i + 3
    status = swPartDoc.InsertCurveFileEnd
    if prev_spline_name == "START":
        feature = swPartDoc.FirstFeature
    else:
        feature = swPartDoc.FeatureByName(prev_spline_name)
    last_feature = None  # Initialize to store the last feature found
    # Traverse through all features until the last one
    while feature:
        if feature.GetTypeName2 in ["CurveInFile", "Curve"]:        
            last_feature = feature  # Update the last feature found
            #print(last_feature.Name)
        feature = feature.GetNextFeature 
    return last_feature.Name

def select_point_in_3d_sketch(prev_feature_name):
    global swSketchManager, swPartDoc
    # Iterate through features to find a 3D Sketch
    feature = swPartDoc.FeatureByName(prev_feature_name)
    last_feature = None
    while feature:
        if feature.GetTypeName2 == "3DProfileFeature":
            last_feature = feature
        feature = feature.GetNextFeature
    sketch = last_feature.GetSpecificFeature2
    # Make sure the sketch is the one with only one point
    sketchPoints = sketch.GetSketchPoints2
    if len(sketchPoints) == 1:  # Only one point in the sketch
        point = sketchPoints[0]
        status = point.SelectByMark(True, 0)
    return status


def create_circular_sweep_profile(first_two_spline_points, spline_name, Sr):
    global swSketchManager, swPartDoc
    P0 = first_two_spline_points[0:3]
    P1 = first_two_spline_points[3:6]
    lx, ly, lz = P1[0] - P0[0], P1[1] - P0[1],P1[2] - P0[2]
    # Normalize the tangent vector
    norm_t = math.sqrt(lx**2 + ly**2 + lz**2)
    lx, ly, lz = lx/norm_t, ly/norm_t, lz/norm_t
    # Dummy vector to get a cross product
    dx, dy, dz = 1, 1, 0
    nx = ly*dz - lz*dy
    ny = lz*dx - lx*dz
    nz = lx*dy - ly*dx
    norm_n = math.sqrt(nx**2 + ny**2 + nz**2)
    nx, ny, nz = nx/norm_n, ny/norm_n, nz/norm_n
    # Create a point at starting location and get its name
    swSketchManager.AddToDB = True
    swSketchManager.DisplayWhenAdded = False
    swSketchManager.Insert3DSketch(True)
    xv, yv, zv = P0[0], P0[1], P0[2]
    swSketchManager.CreatePoint(xv, yv, zv)
    Sketch = swSketchManager.ActiveSketch
    SketchName = Sketch.Name
    swSketchManager.Insert3DSketch(True)
    swSketchManager.AddToDB = False
    swSketchManager.DisplayWhenAdded = True
    Point_select = select_point_in_3d_sketch(spline_name)
    boolstatus = swPartDoc.Extension.SelectByID2(spline_name, "REFERENCECURVES", 0, 0, 0, True, 1, pythoncom.Nothing, 0)
    refplane = swPartDoc.FeatureManager.InsertRefPlane(4, 0, 2, 0, 0, 0)
    status = swPartDoc.Extension.SelectByID2(refplane.Name, "PLANE", 0, 0, 0, False, 0, pythoncom.Nothing, 0)
    swPartDoc.BlankRefGeom()
    status = swPartDoc.Extension.SelectByID2(refplane.Name, "PLANE", 0, 0, 0, False, 0, pythoncom.Nothing, 0)
    swSketchManager.AddToDB = True
    swSketchManager.DisplayWhenAdded = False
    swSketchManager.InsertSketch(True)
    Sketch = swSketchManager.ActiveSketch
    SketchName = Sketch.Name
    circle = swSketchManager.CreateCircleByRadius(0, 0, 0, Sr/1000)
    swPartDoc.ClearSelection2(True)
    swSketchManager.InsertSketch(True)
    swSketchManager.AddToDB = False
    swSketchManager.DisplayWhenAdded = True
    boolstatus = swPartDoc.Extension.SelectByID2(SketchName, "SKETCH", 0, 0, 0, False, 1, pythoncom.Nothing, 0)
    boolstatus = swPartDoc.Extension.SelectByID2(spline_name, "REFERENCECURVES", 0, 0, 0, True, 4, pythoncom.Nothing, 0)
    swFeatData = swPartDoc.FeatureManager.CreateDefinition(17)
    swFeatData.AdvancedSmoothing = False
    swFeatData.AlignWithEndFaces = 0
    swFeatData.AutoSelect = True
    swFeatData.D1ReverseTwistDir = False
    swFeatData.Direction = -1
    swFeatData.EndTangencyType = 0
    swFeatData.FeatureScope = True
    swFeatData.MaintainTangency = False
    swFeatData.Merge = False
    swFeatData.MergeSmoothFaces = True
    swFeatData.PathAlignmentType = 10
    swFeatData.StartTangencyType = 0
    swFeatData.ThinFeature = False
    swFeatData.ThinWallType = 0
    swFeatData.TwistControlType = 0
    swFeatData.SetTwistAngle(0)
    swFeatData.SetWallThickness(True, 0)
    sweep_feature = swPartDoc.FeatureManager.CreateFeature(swFeatData)
    return sweep_feature.Name

def generate_ref_pts(ref_pt):
    # Create a point at starting location and get its name
    global swSketchManager, swPartDoc
    swSketchManager.AddToDB = True
    swSketchManager.DisplayWhenAdded = False
    swSketchManager.Insert3DSketch(True)
    xv, yv, zv = ref_pt[0], ref_pt[1], ref_pt[2]
    swSketchManager.CreatePoint(xv, yv, zv)
    Sketch = swSketchManager.ActiveSketch
    SketchName = Sketch.Name
    swSketchManager.Insert3DSketch(False)
    swSketchManager.AddToDB = False
    swSketchManager.DisplayWhenAdded = True
    return SketchName

def generate_axis_between_pts(ref_pt_1, ref_pt_2):
    # Create a point at starting location and get its name
    global swSketchManager, swPartDoc
    status = swPartDoc.Extension.SelectByID2("Point1@" + ref_pt_1, "EXTSKETCHPOINT", 0, 0, 0, True, 0, pythoncom.Nothing, 0)
    status = swPartDoc.Extension.SelectByID2("Point1@" + ref_pt_2, "EXTSKETCHPOINT", 0, 0, 0, True, 0, pythoncom.Nothing, 0)
    status = swPartDoc.InsertAxis2(True)
    feature = swPartDoc.FeatureByName(ref_pt_2)
    last_feature = None
    while feature:
        if feature.GetTypeName2 in ["RefAxis"]:        
            last_feature = feature  # Update the last feature found
        feature = feature.GetNextFeature 
    return last_feature.Name

def generate_reverse_pattern(body, axis, z_axis, z_shift, theta_shift):
    global swSketchManager, swPartDoc
    swPartDoc.ClearSelection2(True)
    boolstatus = swPartDoc.Extension.SelectByID2(axis, "AXIS", 0, 0, 0, False, 1, pythoncom.Nothing, 0)
    boolstatus = swPartDoc.Extension.SelectByID2(body, "SOLIDBODY", 0, 0, 0, True, 256, pythoncom.Nothing, 0)
    swFeatMgr = swPartDoc.FeatureManager
    swFeatData = swFeatMgr.CreateDefinition(5)
    swFeatData.Direction2 = False
    swFeatData.EqualSpacing = False
    swFeatData.GeometryPattern = False
    swFeatData.ReverseDirection = False
    swFeatData.Spacing = 3.1415926535898
    swFeatData.TotalInstances = 2
    swFeatData.VarySketch = False
    swFeat = swFeatMgr.CreateFeature(swFeatData)
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByID2(swFeat.Name, "SOLIDBODY", 0, 0, 0, False, 1, pythoncom.Nothing, 0)
    translated_feat = swPartDoc.FeatureManager.InsertMoveCopyBody2(0, 0, z_shift, 0, 0, 0, 0, 0, 0, 0, False, 1)    
    swPartDoc.ClearSelection2(True)
    swPartDoc.Extension.SelectByID2(translated_feat.name, "SOLIDBODY", 0, 0, 0, False, 1, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2(z_axis, "AXIS", 0, 0, 0, True, 2, pythoncom.Nothing, 0)
    swPartDoc.Extension.SelectByID2(axis, "AXIS", 0,0,0, True, 4, pythoncom.Nothing, 0)
    rev_feat = swPartDoc.FeatureManager.InsertMoveCopyBody2(0, 0, 0, 0, 0, 0, 0, theta_shift, 0, 0, False, 1)
    return rev_feat.Name


def create_line(first_point,last_point):
    global swSketchManager
    arg1 = win32com.client.VARIANT(pythoncom.VT_R8, first_point[0])
    arg2 = win32com.client.VARIANT(pythoncom.VT_R8, first_point[1])
    arg3_6 = win32com.client.VARIANT(pythoncom.VT_R8, 0)
    arg4 = win32com.client.VARIANT(pythoncom.VT_R8, last_point[0])
    arg5 = win32com.client.VARIANT(pythoncom.VT_R8, last_point[1])
    swSketchManager.CreateLine(arg1,arg2,arg3_6,arg4,arg5,arg3_6)
    return "Created Line"

def end_current_sketch():
    global swSketchManager
    swSketchManager.InsertSketch(True) #Exit the sketch 
    return "sketch exited"

def create_new_sketch():
    global swPartDoc, swSketchManager
    # Select the front plane
    frontPlane = swPartDoc.FeatureByName("Front Plane")
    frontPlane.Select2(False, -1)

def get_active_sketch_name():
    global swSketchManager
    sketch = swSketchManager.ActiveSketch
    sk_feature = sketch
    return sk_feature.Name

def delete_feature(feat_name):
    global swPartDoc, sw
    boolstatus = swPartDoc.Extension.SelectByID2(feat_name, "BODYFEATURE", 0, 0, 0, False, 0, pythoncom.Nothing, 0)
    swPartDoc.EditDelete()
    return True

def start():
    global sw, swPartDoc, swSketchManager
    sw = connectToSW()
    sw.Visible = False
    sw.SetUserPreferenceToggle(706, False)

# Create a new part
    swPartDoc = sw.NewDocument("C:\\ProgramData\\SolidWorks\\SOLIDWORKS 2024\\templates\\Part.prtdot",0,0,0)
    if not swPartDoc:
        return "Failed to create a new part."

# Select the front plane
    frontPlane = swPartDoc.FeatureByName("Front Plane")
    frontPlane.Select2(False, -1)

# Insert a new sketch
    swSketchManager = swPartDoc.SketchManager
    return "started"
    
if __name__ == '__main__':
    startSW()
    sw = connectToSW()
    model = openFile(sw, "path\\to\\SWfile.sldprt")
    updatePrt(model)