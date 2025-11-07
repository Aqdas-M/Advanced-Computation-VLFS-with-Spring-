module MyTime

using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Visualization:paraview_collection
using Gridap.ODEs
# using Gridap.ODEs.TransientFETools
# using Plots
# using Parameters
using WriteVTK

  name::String = "my_time_domain"
  nx::Int = 20
  ny::Int = 5
  order::Int = 2 #4
  ξ::Float64 = 0.0
  vtk_output::Bool = true

  # Fixed parameters
  Lb = 12.5
  mᵨ = 8.36
  EI₁ = 47100.0
  EI₂ = 47100.0
  β = 0.2
  H = 1.1
  α = 0.249
  kₛ = 1e2

  # Domain size
  Ld = Lb # damping zone length
  LΩ = 2Ld + 2Lb
  x₀ = 0.0
  xdᵢₙ = x₀ + Ld
  xb₀ = xdᵢₙ + Lb/2
  xbⱼ = xb₀ + β*Lb
  xb₁ = xb₀ + Lb
  xdₒᵤₜ = LΩ - Ld
  @show Ld
  @show LΩ
  @show x₀
  @show xdᵢₙ
  @show xb₀
  @show xbⱼ
  @show xb₁
  @show xdₒᵤₜ

  # Physics
  g = 9.81
  ρ = 1025
  d₀ = mᵨ/ρ
  a₁ = EI₁/ρ
  a₂ = EI₂/ρ
  # kᵣ = ξ*a₁/Lb

  # wave properties
  λ = α*Lb
  k = 2π/λ
  ω = sqrt(g*k*tanh(k*H))
  T = 2π/ω
  η₀ = 0.01
  ηᵢₙ(x,t) = η₀*cos(k*x[1]-ω*t)
  ϕᵢₙ(x,t) = (η₀*ω/k)*(cosh(k*x[2]) / sinh(k*H))*sin(k*x[1]-ω*t)
  vᵢₙ(x,t) = -(η₀*ω)*(cosh(k*x[2]) / sinh(k*H))*cos(k*x[1]-ω*t)
  vzᵢₙ(x,t) = ω*η₀*sin(k*x[1]-ω*t)
  ηᵢₙ(t::Real) = x -> ηᵢₙ(x,t)
  ϕᵢₙ(t::Real) = x -> ϕᵢₙ(x,t)
  vᵢₙ(t::Real) = x -> vᵢₙ(x,t)
  vzᵢₙ(t::Real) = x -> vzᵢₙ(x,t)

  # Time stepping
  γₜ = 0.5
  βₜ = 0.25
  t₀ = 0.0
  Δt = T/40
  tf = 50*T#/λ_factor
  ∂uₜ_∂u = γₜ/(βₜ*Δt)
  ∂uₜₜ_∂u = 1/(βₜ*Δt^2)

  # Numerics constants
  nx_total = Int(ceil(nx/β)*ceil(LΩ/Lb))
  h = LΩ / nx_total
  γ = 1.0*order*(order-1)/h
  βₕ = 0.5
  αₕ = ∂uₜ_∂u/g * (1-βₕ)/βₕ

  # Damping
  μ₀ = 2.5
  μ₁ᵢₙ(x::VectorValue) = μ₀*(1.0 - sin(π/2*(x[1])/Ld))
  μ₁ₒᵤₜ(x::VectorValue) = μ₀*(1.0 - cos(π/2*(x[1]-xdₒᵤₜ)/Ld))
  μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
  μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k
  ηd(t) = x -> μ₂ᵢₙ(x)*ηᵢₙ(x,t)
  ∇ₙϕd(t) = x -> μ₁ᵢₙ(x)*vzᵢₙ(x,t)

  # Fluid model
  domain = (x₀, LΩ, 0.0, H)
  partition = (nx_total,ny)
  function f_y(x)
    if x == H
        return H
    end
    i = x / (H/ny)
    return H-H/(2.5^i)
  end
  map(x) = VectorValue(x[1], f_y(x[2]))
  𝒯_Ω = CartesianDiscreteModel(domain,partition,map=map)

  # Labelling
  labels_Ω = get_face_labeling(𝒯_Ω)
  add_tag_from_tags!(labels_Ω,"surface",[3,4,6])   # assign the label "surface" to the entity 3,4 and 6 (top corners and top side)
  add_tag_from_tags!(labels_Ω,"bottom",[1,2,5])    # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
  add_tag_from_tags!(labels_Ω,"inlet",[7])         # assign the label "inlet" to the entity 7 (left side)
  add_tag_from_tags!(labels_Ω,"outlet",[8])        # assign the label "outlet" to the entity 8 (right side)
  add_tag_from_tags!(labels_Ω, "water", [9])       # assign the label "water" to the entity 9 (interior)
  # Triangulations
  Ω = Interior(𝒯_Ω)
  Γ = Boundary(𝒯_Ω,tags="surface")
  Γin = Boundary(𝒯_Ω,tags="inlet")

  # Auxiliary functions
  function is_beam1(xs) # Check if an element is inside the beam1
    n = length(xs)
    x = (1/n)*sum(xs)
    (xb₀ <= x[1] <= xbⱼ ) * ( x[2] ≈ H)
  end
  function is_beam2(xs) # Check if an element is inside the beam2
    n = length(xs)
    x = (1/n)*sum(xs)
    (xbⱼ <= x[1] <= xb₁ ) * ( x[2] ≈ H)
  end
  function is_damping1(xs) # Check if an element is inside the damping zone 1
    n = length(xs)
    x = (1/n)*sum(xs)
    (x₀ <= x[1] <= xdᵢₙ ) * ( x[2] ≈ H)
  end
  function is_damping2(xs) # Check if an element is inside the damping zone 2
    n = length(xs)
    x = (1/n)*sum(xs)
    (xdₒᵤₜ <= x[1] ) * ( x[2] ≈ H)
  end
  function is_beam_boundary(xs) # Check if an element is on the beam boundary
    is_on_xb₀ = [x[1]≈xb₀ for x in xs] # array of booleans of size the number of points in an element (for points, it will be an array of size 1)
    is_on_xb₁ = [x[1]≈xb₁ for x in xs]
    element_on_xb₀ = minimum(is_on_xb₀) # Boolean with "true" if at least one entry is true, "false" otherwise.
    element_on_xb₁ = minimum(is_on_xb₁)
    element_on_xb₀ | element_on_xb₁ # Return "true" if any of the two cases is true
  end
  function is_a_joint(xs) # Check if an element is a joint
    is_on_xbⱼ = [x[1]≈xbⱼ && x[2]≈H for x in xs] # array of booleans of size the number of points in an element (for points, it will be an array of size 1)
    element_on_xbⱼ = minimum(is_on_xbⱼ) # Boolean with "true" if at least one entry is true, "false" otherwise.
    element_on_xbⱼ
  end

  # Beam triangulations
  xΓ = get_cell_coordinates(Γ)
  Γb1_to_Γ_mask = lazy_map(is_beam1,xΓ)
  Γb2_to_Γ_mask = lazy_map(is_beam2,xΓ)
  Γd1_to_Γ_mask = lazy_map(is_damping1,xΓ)
  Γd2_to_Γ_mask = lazy_map(is_damping2,xΓ)
  Γb1_to_Γ = findall(Γb1_to_Γ_mask)
  Γb2_to_Γ = findall(Γb2_to_Γ_mask)
  Γd1_to_Γ = findall(Γd1_to_Γ_mask)
  Γd2_to_Γ = findall(Γd2_to_Γ_mask)
  Γf_to_Γ = findall(!,Γb1_to_Γ_mask .| Γb2_to_Γ_mask .| Γd1_to_Γ_mask .| Γd2_to_Γ_mask)
  Γη_to_Γ = findall(Γb1_to_Γ_mask .| Γb2_to_Γ_mask )
  Γκ_to_Γ = findall(!,Γb1_to_Γ_mask .| Γb2_to_Γ_mask )
  Γb1 = Triangulation(Γ,Γb1_to_Γ)
  Γb2 = Triangulation(Γ,Γb2_to_Γ)
  Γd1 = Triangulation(Γ,Γd1_to_Γ)
  Γd2 = Triangulation(Γ,Γd2_to_Γ)
  Γfs = Triangulation(Γ,Γf_to_Γ)
  Γη = Triangulation(Γ,Γη_to_Γ)
  Γκ = Triangulation(Γ,Γκ_to_Γ)
  Λb1 = Skeleton(Γb1)
  Λb2 = Skeleton(Γb2)

  # Construct the mask for the joint
  Γ_mask_in_Ω_dim_0 = get_face_mask(labels_Ω,"surface",0)
  grid_dim_0_Γ = GridPortion(Grid(ReferenceFE{0},𝒯_Ω),Γ_mask_in_Ω_dim_0)
  xΓ_dim_0 = get_cell_coordinates(grid_dim_0_Γ)
  Λj_to_Γ_mask = lazy_map(is_a_joint,xΓ_dim_0)
  Λj = Skeleton(Γ,Λj_to_Γ_mask)


  if vtk_output == true
    filename = "Output/time/k/results/"*name
    writevtk(Ω,filename*"_O")
    writevtk(Γ,filename*"_G")
    writevtk(Γb1,filename*"_Gb1")
    writevtk(Γb2,filename*"_Gb2")
    writevtk(Γd1,filename*"_Gd1")
    writevtk(Γd2,filename*"_Gd2")
    writevtk(Γfs,filename*"_Gfs")
    writevtk(Λb1,filename*"_L1")
    writevtk(Λb2,filename*"_L2")
    writevtk(Λj,filename*"_Lj")
  end

  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓb1 = Measure(Γb1,degree)
  dΓb2 = Measure(Γb2,degree)
  dΓd1 = Measure(Γd1,degree)
  dΓd2 = Measure(Γd2,degree)
  dΓfs = Measure(Γfs,degree)
  dΓin = Measure(Γin,degree)
  dΛb1 = Measure(Λb1,degree)
  dΛb2 = Measure(Λb2,degree)
  dΛj = Measure(Λj,degree)

  # Normals
  nΛb1 = get_normal_vector(Λb1)
  nΛb2 = get_normal_vector(Λb2)
  nΛj = get_normal_vector(Λj)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_Ω = TestFESpace(Ω, reffe, conformity=:H1)
  V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1)
  V_Γη = TestFESpace(Γη, reffe, conformity=:H1)
  U_Ω = TransientTrialFESpace(V_Ω)
  U_Γκ = TransientTrialFESpace(V_Γκ)
  U_Γη = TransientTrialFESpace(V_Γη)
  X = TransientMultiFieldFESpace([U_Ω,U_Γκ,U_Γη])
  Y = MultiFieldFESpace([V_Ω,V_Γκ,V_Γη])

  # Weak form
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  m(t,(ϕₜₜ,κₜₜ,ηₜₜ),(w,u,v)) = ∫( d₀*ηₜₜ*v )dΓb1 + ∫( d₀*ηₜₜ*v )dΓb2
  c(t,(ϕₜ,κₜ,ηₜ),(w,u,v)) = ∫( βₕ*ϕₜ*(u + αₕ*w) - κₜ*w )dΓfs +
                        ∫( βₕ*ϕₜ*(u + αₕ*w) - κₜ*w )dΓd1 +
                        ∫( βₕ*ϕₜ*(u + αₕ*w) - κₜ*w )dΓd2 +
                        ∫( ϕₜ*v - ηₜ*w )dΓb1 +
                        ∫( ϕₜ*v - ηₜ*w )dΓb2
  a(t,(ϕ,κ,η),(w,u,v)) =  ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
                        ∫(  βₕ*(u + αₕ*w)*(g*κ) )dΓfs   +
                        ∫(  βₕ*(u + αₕ*w)*(g*κ) - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1     +
                        ∫(  βₕ*(u + αₕ*w)*(g*κ) - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2    +
                        ∫(  ( v*(g*η) + a₁*Δ(v)*Δ(η) ) )dΓb1  +
                        ∫(  ( v*(g*η) + a₂*Δ(v)*Δ(η) ) )dΓb2  +
                        ∫(  a₁ * ( - jump(∇(v)⋅nΛb1) * mean(Δ(η)) - mean(Δ(v)) * jump(∇(η)⋅nΛb1) + γ*( jump(∇(v)⋅nΛb1) * jump(∇(η)⋅nΛb1) ) ) )dΛb1 +
                        ∫(  a₂ * ( - jump(∇(v)⋅nΛb2) * mean(Δ(η)) - mean(Δ(v)) * jump(∇(η)⋅nΛb2) + γ*( jump(∇(v)⋅nΛb2) * jump(∇(η)⋅nΛb2) ) ) )dΛb2 +
                        ∫(  a₁ * ( - jump(∇(v)⋅nΛj) * mean(Δ(η)) - mean(Δ(v)) * jump(∇(η)⋅nΛj) + γ*( jump(∇(v)⋅nΛj) * jump(∇(η)⋅nΛj) ) ) )dΛj +
                        # ∫(  (jump(∇(v)⋅nΛj) * kᵣ * jump(∇(η)⋅nΛj)) )dΛj
                        ∫(  ( kₛ * mean(η) * mean(v)) )dΛj
  b(t,(w,u,v)) =  ∫( w*vᵢₙ(t) )dΓin - ∫( ηd(t)*w - ∇ₙϕd(t)*(u + αₕ*w) )dΓd1

  # ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  # m((ϕₜₜ,κₜₜ,ηₜₜ),(w,u,v)) = ∫( d₀⋅ηₜₜ⋅v )dΓb1 + ∫( d₀⋅ηₜₜ⋅v )dΓb2
  # c((ϕₜ,κₜ,ηₜ),(w,u,v)) = ∫( βₕ⋅ϕₜ⋅(u + αₕ⋅w) - κₜ⋅w )dΓfs +
  #                       ∫( βₕ⋅ϕₜ⋅(u + αₕ⋅w) - κₜ⋅w )dΓd1 +
  #                       ∫( βₕ⋅ϕₜ⋅(u + αₕ⋅w) - κₜ⋅w )dΓd2 +
  #                       ∫( ϕₜ⋅v - ηₜ⋅w )dΓb1 +
  #                       ∫( ϕₜ⋅v - ηₜ⋅w )dΓb2
  # a((ϕ,κ,η),(w,u,v)) =  ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
  #                       ∫(  βₕ⋅(u + αₕ⋅w)⋅(g⋅κ) )dΓfs   +
  #                       ∫(  βₕ⋅(u + αₕ⋅w)⋅(g⋅κ) - μ₂ᵢₙ⋅κ⋅w + μ₁ᵢₙ⋅∇ₙ(ϕ)⋅(u + αₕ⋅w) )dΓd1     +
  #                       ∫(  βₕ⋅(u + αₕ⋅w)⋅(g⋅κ) - μ₂ₒᵤₜ⋅κ⋅w + μ₁ₒᵤₜ⋅∇ₙ(ϕ)⋅(u + αₕ⋅w) )dΓd2    +
  #                       ∫(  ( v⋅(g⋅η) + a₁⋅Δ(v)⋅Δ(η) ) )dΓb1  +
  #                       ∫(  ( v⋅(g⋅η) + a₂⋅Δ(v)⋅Δ(η) ) )dΓb2  +
  #                       ∫(  a₁⋅( - jump(∇(v)⋅nΛb1)⋅mean(Δ(η)) - mean(Δ(v))⋅jump(∇(η)⋅nΛb1) + γ⋅( jump(∇(v)⋅nΛb1)⋅jump(∇(η)⋅nΛb1) ) ) )dΛb1 +
  #                       ∫(  a₂⋅( - jump(∇(v)⋅nΛb2)⋅mean(Δ(η)) - mean(Δ(v))⋅jump(∇(η)⋅nΛb2) + γ⋅( jump(∇(v)⋅nΛb2)⋅jump(∇(η)⋅nΛb2) ) ) )dΛb2 +
  #                       # ∫(  (jump(∇(v)⋅nΛj)⋅kᵣ⋅jump(∇(η)⋅nΛj)) )dΛj
  #                       ∫(  ( kₛ⋅mean(η)⋅mean(v)) )dΛj
  # b(t,(w,u,v)) =  ∫( w⋅vᵢₙ(t) )dΓin - ∫( ηd(t)⋅w - ∇ₙϕd(t)⋅(u + αₕ⋅w) )dΓd1

  # Solution
  # a c m 
  # op = TransientConstantMatrixFEOperator(m,c,a,b,X,Y)
  op = TransientLinearFEOperator((a,c,m),b,X,Y)
  ls = LUSolver()
  ode_solver = Newmark(ls,Δt,γₜ,βₜ)

  # Initial solution
  x₀ = interpolate_everywhere([0.0,0.0,0.0],X(0.0))
  v₀ = interpolate_everywhere([0.0,0.0,0.0],X(0.0))
  a₀ = interpolate_everywhere([0.0,0.0,0.0],X(0.0))

  xₜ = solve(ode_solver,op,t₀,tf,(x₀,v₀,a₀))

  if vtk_output == true
    pvd_Ω = paraview_collection(filename * "_O_solution", append=false)
    pvd_Γκ = paraview_collection(filename * "_Gk_solution", append=false)
    pvd_Γη = paraview_collection(filename * "_Ge_solution", append=false)
  end

  # Postprocess
  xy_cp = get_cell_points(get_fe_dof_basis(V_Γη)).cell_phys_point
  x_cp = [[xy_ij[1] for xy_ij in xy_i] for xy_i in xy_cp]
  p = sortperm(x_cp[1])
  x_cp_sorted = [x_i[p] for x_i in x_cp]
  xs = [(x_i-1.5*Lb)/Lb for x_i in vcat(x_cp_sorted...)]


  ts = Float64[]
  ηxps_t = []
  for (tₙ,(ϕₕ,κₕ,ηₕ)) in xₜ
    println("t = $tₙ")
    η_cdv = get_cell_dof_values(ηₕ)
    η_cdv_sorted = [η_i[p] for η_i in η_cdv]
    η_rel_xs = [abs(η_i)/η₀ for η_i in vcat(η_cdv_sorted...)]
    push!(ηxps_t,η_rel_xs)
    push!(ts,tₙ)

    if vtk_output == true
      pvd_Ω[tₙ] = createvtk(Ω,filename * "_O_solution" * "_$tₙ.vtu",cellfields = ["phi" => ϕₕ])#,nsubcells=10)
      pvd_Γκ[tₙ] = createvtk(Γκ,filename * "_Gk_solution" * "_$tₙ.vtu",cellfields = ["kappa" => κₕ],nsubcells=10)
      pvd_Γη[tₙ] = createvtk(Γη,filename * "_Ge_solution" * "_$tₙ.vtu",cellfields = ["eta" => ηₕ],nsubcells=10)
    end
  end

  if vtk_output == true
    vtk_save(pvd_Ω)
    vtk_save(pvd_Γκ)
    vtk_save(pvd_Γη)
  end

  return (ts,xs,ηxps_t)
end
