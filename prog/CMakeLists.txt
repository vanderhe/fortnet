add_subdirectory(fortnet)

# Propagate PKG_CONFIG_LIBS upwards
set(PKG_CONFIG_LIBS "${PKG_CONFIG_LIBS}" PARENT_SCOPE)
