:orphan:

More information about containers
=================================

The ``Module`` object is defined as a named ``Parameters`` container.
Its naming is stored within its ``mod_name`` attribute.
For instance, the example shown `in this section </cards-python>`_ may be rewritten as:

.. code-block:: python

   import Config.Core as cepgen
   module = cepgen.Parameters(mod_name = 'my_first_module', foo = 'bar')

or, using standard Python containers:

.. code:: python

   module = {
       'mod_name': 'my_first_module',
       'foo': 'bar'
   }

Among the features introduced by these two objects, one may quote the capability of deep-cloning one while changing one or several attributes in one pass.
For instance, using the former definition of ``module``:

.. code:: python

   module_bis = module.clone('my_cloned_module',
       foo2 = 42,
       foo3 = cepgen.Parameters(
           nestedFoo = 'nestedBar',
       )
   )

is equivalent to defining a ``module_bis`` container as:

.. code:: python

   module_bis = cepgen.Module('my_cloned_module',
       foo = 'bar',
       foo2 = 42,
       foo3 = cepgen.Parameters(
           nestedFoo = 'nestedBar',
       ),
   )

