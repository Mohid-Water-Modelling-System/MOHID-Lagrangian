module tracer 
    ! Wrapper to hold all libtracer modules 
    
    use bspline_module, only : bspline_3d 
    
    use tracer_precision
    use tracer_interp 
    
    use tracer3D 
    use tracer2D 

    use tracer_io 

end module tracer 
