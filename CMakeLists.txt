###################################################################
# Package: SuperNtuple
###################################################################

# declare the package name
atlas_subdir( SuperNtuple )

# extra deps based on the build environment
set( extra_deps )
set( extra_libs )

# declare the SuperNtuple package dependencies
atlas_depends_on_subdirs(
    PUBLIC
    Superflow
    DileptonTriggerTool
    SusyNtuple
    hhTruthWeightTools
    HHLambdaTool
    ${extra_deps}
)

# extra dependencies
find_package( ROOT COMPONENTS Gpad Graf Graf3d Core Tree MathCore Hist RIO )
find_package( AnalysisBaseExternals REQUIRED )
find_package( lwtnn REQUIRED )

# libraries
atlas_add_library( SuperNtupleLib
    SuperNtuple/*.h Root/*.cxx
    PUBLIC_HEADERS SuperNtuple
    PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} ${LWTNN_INCLUDE_DIRS}
    LINK_LIBRARIES SuperflowLib SusyNtupleLib DileptonTriggerToolLib
        hhTruthWeightToolsLib HHLambdaToolLib 
        ${ROOT_LIBRARIES} ${LWTNN_LIBRARIES}
)

# executables
set( extra_libs )

# define a function to call the add_executable macro on our executables
function(SuperNtupleExec filename)
    set(execname)
    get_filename_component(execname ${filename} NAME_WE)
    atlas_add_executable( ${execname} "util/${execname}.cxx"
        INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} ${LWTNN_INCLUDE_DIRS}
        LINK_LIBRARIES ${ROOT_LIBRARIES} ${LWTNN_LIBRARIES} SuperNtupleLib ${extra_libs}
    )
endfunction(SuperNtupleExec)

# build everything in util/ as an executable
file(GLOB files "util/*.cxx")
foreach(file ${files})
    SuperNtupleExec(${file})
endforeach()

atlas_install_data( data/* )
