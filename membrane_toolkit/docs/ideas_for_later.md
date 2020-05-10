Key Classes
-----------

The primary unit of analysis for transport calcs is a differential membrane
slice comprising 5 objects (probably 3 classes):

Bulk | Interface | Membrane | Interface | Bulk

Each of these "phase" classes (bulk, interface, Membrane) will be a base class
that can be customized e.g. IXMembrane, ROMembrane, AqueousSolution, OrganicSolution,
etc. Bulk could even be a vapor class

All of these classes need to have a consistent interface to access 1) composition
2) state variables like pressure, density, etc., and 3) (electro) chemical potential.

the customization will specify e.g. how to calculate chemical potential from composition
AqueousSolution might use D-H or Pitzer models for activity and density, whereas a
Membrane class might use Manning theory. 

The member methods (e.g. get_activity_coefficient_manning) should be usable by themselves

How to handle driving forces? I guess you'll specify these at the boundaries (bulk)
and allow the intermediate phases to be solved?