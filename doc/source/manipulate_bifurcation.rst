.. title:: Tutorial: Manipulate bifurcation

================================
Tutorial: Manipulate bifurcation
================================

The script ``move_branches.py`` performes an objective rotation of two daughter branches, moving them a given angle along the bifurcation plane. 
In this implementation, positive and negative angles rotate the branches upward and downward, respectively.
The implementation rotates both branches by default, but the user can specificy that only rotate a single branch. 
This is achieved by setting the argument ``keep-fixed1`` or ``keep-fixed2`` to **True**, to leave the first or second branch, respectively.

To perform rotation of the daughter branches, run the following generalized command::
    
    python move_branches.py --dir_path [PATH_TO_CASES] --case [CASENAME] --angle [ANGLE] --lower True


.. figure:: angle_updown.png

  Figure 5: Rotation of daughter branches, in both a widening and narrowing of the bifurcation angle. 

In addition, the user includes two options for reconstruction of the bifurcation, determined by the parameters ``bif`` and ``lower`` set to **True**.
The ``lower`` parameter creates a more realistic looking bifurcation, while the ``bif`` parameter creates a straight segment between the daughter branches.
[\\]: <> (Cite or refer to reconstruction paper)
A comparison is shown below, where the straight segment is created by running::

    python move_branches.py --dir_path [PATH_TO_CASES] --case [CASENAME] --angle [ANGLE] --bif True

.. figure:: angle_bif.png

  Figure 6: Rotation of daughter branches with different reconstruction of the bifurcation.


