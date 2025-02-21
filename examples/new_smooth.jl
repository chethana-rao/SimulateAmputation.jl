using Comodo
using GLMakie
using GeometryBasics
using FileIO
using Statistics
using LinearAlgebra
using Geogram

# 
Z_cut_offset_top = 200.0
Z_cut_level_skin = 100.0 #(100-160)
Z_cut_level_tibia = Z_cut_level_skin+50
Z_cut_level_fibula = Z_cut_level_skin+10

Z_thickness_distal = 20 #distal thickness can be changed

pointSpacing = 6.0 
fileName_set = (
     "/home/simulimb/Che_simamp/SimulateAmputation.jl/Che_assets/FMA7163_right_leg_isolated_remesh_25k.stl",
"/home/simulimb/Che_simamp/SimulateAmputation.jl/Che_assets/FMA24474.stl",
"/home/simulimb/Che_simamp/SimulateAmputation.jl/Che_assets/FMA24486.stl",
"/home/simulimb/Che_simamp/SimulateAmputation.jl/Che_assets/FMA24477.stl",
"/home/simulimb/Che_simamp/SimulateAmputation.jl/Che_assets/FMA24480.stl")

nameSet = ("skin","femur","patella","tibia","fibula")
indSkin = 1
indFemur = 2
indPatella = 3
indTibia = 4
indFibula = 5

nameSetCut_top = ("skin","femur")
# nameSetCut_bottom = ("skin","femur","tibia","fibula")

n = normalizevector(Vec{3, Float64}(0.0,0.0,1.0))# Cutting plane normal
snapTolerance = 1e-6

# FMAID 	English name
# 7163  	Skin
# 16586 	Right hip bone
# 24474 	Right femur
# 24486     Right patella
# 24477 	Right tibia
# 24480 	Right fibula

FF = []
VV = []
for i in 1:1:length(fileName_set)
    M = load(fileName_set[i])
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    V = [Point{3,Float64}(v) for v in V]
    
    F,V,_,_ = mergevertices(F,V; pointSpacing=pointSpacing)

    # Remesh evenly using desired point spacing

    # Determine number of points to request to closely match desired point spacing
    np = spacing2numvertices(F,V,pointSpacing) 

    # Remesh with desired number of points 
    F,V = ggremesh(F,V; nb_pts=np,suppress = true, cleanup = false)

    push!(FF,F)
    push!(VV,V)
end

FF_ori = deepcopy(FF)
VV_ori = deepcopy(VV)

# Cutting surfaces 
#finding the midpoint of patella 
v_mid_patella = mean(VV[indPatella])

