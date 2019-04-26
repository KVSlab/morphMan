.. title:: Tutorial: Manipulate branch

.. _manipulate_branch:

===========================
Tutorial: Manipulate branch
===========================

The goal of ``manipulate_branch.py`` is to manipulate a specific branch of a
vascular model, either through pure tranlation and rotation, or by completely removing it, as depicted in Figure 1.
In particular, we have defined a branch as any branch branching out of the longest tubular structure,
defined by the physical length of the tube.
The manipulation can be achieved by running ``morphman-branch`` in the terminal, followed by the
respective command line arguments. Alternatively, you can execute the Python script directly,
located in the ``morphman`` subfolder, by typing ``python manipulate_branch.py``. We have also created a
demo folder where we show how to run this tutorial from a Python script, please checkout the code from GitHub to
run the demos.

.. figure:: Branch.png

  Figure 1: An illustration of the desired output from the method.

In this tutorial, we are using the model with
`ID C0002 <http://ecm2.mathcs.emory.edu/aneuriskdata/download/C0002/C0002_models.tar.gz>`_
from the Aneurisk database. For the commands below we assume that there is a
file `./C0002/surface/model.vtp`, relative to where you execute the command.

When using ``morphman-branch``, there are four main settings which can be provided by the user:

 * ``branch_number``: With branches ordered from 1 to N from upstream to downstream relative to the inlet, this number determines which branch is to be manipulated.
 * ``branch_loc``: The point on / closes to the surface, where the selected branch will be placed.
 * ``angle``: How many degrees the manipulated branch will be rotated around the new surface normal vector.
 * ``remove_branch``: Either `True` or `False`, which determines whether to remove the selected branch or not.

Shown in Figure 2 is the result of moving the opthalmic artery of the vascular model to another part of the surface.

.. figure:: manipulate_opthalmic.png

  Figure 2: Translation and rotation of the opthamlic artery,
  causing it to appear elsewhere on the vasular surface model.

To reproduce the surface model where the opthalmic artey has been moved, as shown in Figure 2, run::

    morphman-branch --ifile C0002/surface/model.vtp --ofile C0002/surface/moved_branch.vtp --branch-number 1 --branch-location 21.7 18.1 25.9
  --poly-ball-size 250 250 250

As explained earlier, setting the `branch-number` equal to 1 corresponds to the first branch out of the main tube,
in this case the opthalmic artery of the ICA model.

Shown in Figure 2 is the result of moving the opthalmic artery of the vascular model to another part of the surface,
including a second rotation around the new surface normal vector, which ranges of the angle :math:`\theta \in [0, 2 \pi ]`.

.. figure:: manipulate_opthalmic_rotation.png

  Figure 3: Translation and rotation of the opthamlic artery, followed by rotation around the new surface normal vector.

To reproduce the surface model where the opthalmic artey has been moved, as shown in Figure 3, run::

    morphman-branch --ifile C0002/surface/model.vtp --ofile C0002/surface/moved_and_rotated_branch.vtp --angle 180 --branch-number 1 --branch-location 21.7 18.1 25.9  --poly-ball-size 250 250 250

Notice how the `angle` setting is given in degrees, although converted to radians in the main algorithm.

Finally, in Figure 4 is the result of removing an arbitrary branch.

.. figure:: manipulate_opthalmic_rotation.png

  Figure 4: Removal of an arbitrary branch, in this case the third branch away from the inlet.

To reproduce the surface model where a branch has been removed, as shown in Figure 4, run::

    morphman-branch --ifile C0002/surface/model.vtp --ofile C0002/surface/removed_branch.vtp --remove-branch True --branch-number 3 --poly-ball-size 250 250 250

For additional information, beyond this tutorial, on the script and
input parameters, please run ``morphman-branch -h`` or confer with
the :meth:`curvature_branch.curvature_branch`.
