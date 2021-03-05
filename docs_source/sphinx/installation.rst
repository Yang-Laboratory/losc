============
Installation
============

.. note:: Currently, the library can only be used in linux platform.

----------------------
Dependencies and tools
----------------------

- To build the C++/C parts of the LOSC library, following tools and dependencies
  are necessary:

    - C++ and C compilers (C++11 compliant)
    - `Eigen3 library <https://eigen.tuxfamily.org/dox/>`_
    - `CMake <http://www.cmake.org/download/>`_
    - `Openmp <https://www.openmp.org/>`_
    - System utilities: GNU make, GNU install,

- For the Python parts of the LOSC library, namely the
  :ref:`py_losc <py_losc:LOSC Python library>` and
  :ref:`psi4_losc <psi4_losc:use LOSC library in psi4>` modules,
  there is no need to compile. However, following dependencies are needed
  at the running time:

    - `Python interpreter (3.7) <https://www.python.org/>`_
    - `Numpy <http://www.numpy.org/>`_
    - `The psi4 package <https://psicode.org/>`_ for using the
      ``psi4_losc`` module.

- There is no need to build the documentation by yourself. You can find the
  `online documentation <https://miocbb.github.io/losc>`_. However, if you
  insist to generate the documentation that is corresponding to your local
  version of the source code, you need the following dependencies and tools:

    - `sphinx <https://www.sphinx-doc.org/en/master/#>`_: python documentation
      generator.
    - `doxygen <https://www.doxygen.nl/index.html>`_: C/C++ documentation
      generator.

----------------------------
Compile from the source code
----------------------------

#. you can obtain the source code from `the github repository
   <https://github.com/Miocbb/losc>`_ via ``git clone``:

   .. code-block:: bash

        >>> git clone https://github.com/Miocbb/losc.git

   Doing so requires you to have ``git`` installed in your system. If you do
   not have access to ``git``, you can download the archive of the source
   code from the `the github repository <https://github.com/Miocbb/losc>`_,
   and extract the source manually.

#. Navigate to the top level directory of the source code
   ``{top-level-losc-dir}``. It is suggested to separate the building files
   to the source files with creating a build directory ``{build_dir}``.
   Then run ``cmake`` to configure the building process.

   .. code-block:: bash

        >>> cd {top-level-losc-dir}
        >>> mkdir {build_dir}
        >>> cd {build_dir}
        >>> cmake {top-level-losc-dir}

   You can configure ``cmake`` with following options:

       - ``DCMAKE_C_COMPILER``: The executable path for C compiler
       - ``DCMAKE_CXX_COMPILER``: The executable path for C++ compiler
       - ``DCMAKE_INSTALL_PREFIX``: The prefix path for installing the library.
         Default to ``/usr/local``.
       - ``DOPENMP``: Enable ``openmp`` for parallel threading or not.
         Default to on.
       - ``DBUILD_TEST``: Build the test or not. Default to off.
       - ``DBUILD_DOC``: Build the documentation or not. Default to off.

#. Invoke ``make`` to build and install.

   .. code-block:: bash

      >>> make
      >>> make install # installing is optional.

   or invoke``make`` in parallel with multiple processors.

   .. code-block:: bash

      >>> make -j {number-of-processors}

   **If all you need is the C/C++ parts of the LOSC library, you are okay to
   stop here. If you want to use the Python interface or the LOSC library
   in psi4, please continue.**

#. Configure running time for Python modules (``py_losc`` and ``psi4_losc``).

   Append the path of ``py_losc`` and ``psi4_losc`` module to the ``PYTHONPATH``
   environment variable in your ``~/.zshrc`` or ``~/.bashrc`` file. This is to
   enable the Python interpreter to locate ``py_losc`` and ``psi4_losc`` modules
   and import them successfully at running time.

   If you installed the LOSC library:

   .. code-block:: bash

      export PYTHONPATH=${PYTHONPATH}:{DCMAKE_INSTALL_PREFIX}/liblosc

   If you only build the LOSC library and not install it:

   .. code-block:: bash

      export PYTHONPATH=${PYTHONPATH}:{build_dir}/src

#. Running tests to verify the compilation/installing is optional.

    - To run tests for ``losc`` C++ library, remember to build tests for
      ``losc`` first with ``DBUILD_TEST=On``. Then run the executable losc test
      file.

      .. code-block:: bash

         >>> {build_dir}/tests/losc/losc_test

    - To run tests for ``psi4_losc`` Python module, do the following.

      .. code-block:: bash

         >>> cd {top-level-losc-dir}/tests/psi4_losc
         >>> python3 -m unittest test_scf_losc.py  --verbose

    - There are no tests for ``py_losc`` Python module.

--------------------------
Uninstall the LOSC library
--------------------------

To uninstall the LOSC library, remove the whole installed directory of
LOSC.

.. code-block:: bash

   >>> rm {DCMAKE_INSTALL_PREFIX}/liblosc
