project(
    'quspin-core',
    'cpp',
    version: '0.1.0',
    default_options: ['warning_level=3', 'cpp_std=c++20'],
)

add_project_arguments('-DVERSION_INFO=' + meson.project_version() + '', language: 'cpp')

# libquspin subproject
libquspin = subproject('libquspin')
libquspin_dep = libquspin.get_variable('libquspin_dep')

# python bindings
py = import('python').find_installation(pure: false)
pybind11_dep = dependency('pybind11')

modules = [
    '_version',
    'array',
]

foreach module : modules
    py.extension_module(module,
        'src/quspin_core/' + module + '.cpp',
        install: true,
        subdir: 'quspin_core',
        dependencies : [pybind11_dep, libquspin_dep],
    )
endforeach

python_sources = files('src/quspin_core/__init__.py')


py.install_sources(
  python_sources,
  subdir: 'quspin_core'
)
