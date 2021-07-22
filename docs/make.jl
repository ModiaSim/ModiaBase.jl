using Documenter, ModiaBase

makedocs(
  #modules  = [ModiaBase],
  sitename = "ModiaBase",
  authors  = "Hilding Elmqvist (Mogram) and Martin Otter (DLR-SR)",
  format = Documenter.HTML(prettyurls = false),
  pages    = [
     "Home"   => "index.md",
	 "Tutorial" => "Tutorial.md",
     "Data Structures"        => "OrderedCollections.md",
     "Equation Sorting"       => "EquationSorting.md",
     "Equation Reduction"     => "EquationReduction.md",
     "Transformation to ODE System" => "TransformationToODEs.md",
  ]
)
