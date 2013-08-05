ADD_LIBRARY (RootIoWrapper SHARED event_reader.cpp event_writer.cpp mydict.cxx)
target_link_libraries (RootIoWrapper ${ROOT_LIBRARIES} ${GSL_LIBRARIES})

ADD_CUSTOM_COMMAND(
      OUTPUT mydict.cxx
      COMMAND rootcint
      ARGS -f mydict.cxx -c vector ../../include/ROOT_Ttree_io/myLinkDef.h
)
