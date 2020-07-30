def plot_it(dataArray, output_filename, title, units, plot_type, polygon_file):
    """


    Args:
        dataArray (TYPE): DESCRIPTION.
        output_filename (TYPE): DESCRIPTION.
        plot_type (TYPE, optional): DESCRIPTION. Defaults to 'contourf'.

    Returns:
        None.

    """
    if polygon_file:
        geoDataFrame = gpd.read_file(polygon_file).to_crs({"init": 'EPSG:4326'})
        geoDataFrame['index'] = np.arange(0, geoDataFrame.shape[0])

    time_plot_flag = hastime(dataArray)

    if time_plot_flag:
        time_key = dataArray.dims[0]
        time_size = dataArray.shape[0]
        if time_size < 4:
            nrows = 1
            ncols = time_size
        elif time_size == 4:
            nrows = 2
            ncols = 2
        elif time_size == 6:
            nrows = 2
            ncols = 3
        elif time_size == 8:
            nrows = 2
            ncols = 4
        elif time_size == 12:
            nrows = 3
            ncols = 4
        else:
            time_slice = slice(0, -1, int(time_size/4))
            dataArray = dataArray.isel({time_key: time_slice})
            nrows = 2
            ncols = 2
    else:
        nrows = ncols = 1

    if nrows == ncols:
        figsize = (17, 15)
    elif ncols > nrows:
        figsize = (20, 15)

    fig, axarr = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize,
                              subplot_kw={'projection': ccrs.PlateCarree()})

    # generalization -> make one subplot iterable.
    if not isinstance(axarr, np.ndarray):
        axarr = np.array([axarr])

    vmin, vmax = get_color_lims(dataArray, robust=True)
    #scale_bar_lenght = get_horizontal_scale(dataArray)
    cmap_key = get_cmap_key(vmin, vmax)
    cmap_norm = get_cmap_norm(vmin, vmax)
    time_step = 0

    for ax in axarr.flat:

        if time_plot_flag:
            dataArray_step = dataArray.isel({time_key: time_step})
        else:
            dataArray_step = dataArray

        if polygon_file:
            geoDataFrame = geoDataFrame.merge(dataArray.to_dataframe(),on='index')
            extent = get_polygon_extent(geoDataFrame)
            varName = dataArray.name
        else: 
            extent = get_extent(dataArray)

        #try:
        setup_plot = {
                      'cmap': cmap_key,
                      'ax': ax,
                      'vmin': vmin,
                      'vmax': vmax,
                      'zorder': 1}
                      #'add_colorbar': False}

        if plot_type == 'hist':
             setup_plot ={'ax': ax,
                          'zorder': 1}

        if polygon_file:
            geoDataFrame.plot(column=varName,**setup_plot)
            print('PLOTTING THE GEODATAFRAME')
        else:
            p = getattr(dataArray_step.plot, plot_type)(**setup_plot)

        ax = get_background_map(ax, extent)
        gl = ax.gridlines(draw_labels=True, color='gray', linestyle='--')
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        if np.abs((extent[3]-extent[2]) > ((extent[1]-extent[0]))):
            gl.xlabel_style = {'rotation': 45}

        #scale_bar(ax, ccrs.PlateCarree(), scale_bar_lenght)
        time_step += 1

     #  except Exception as exc:
     #       print('-> WARNING: Could not plot', title)
     #       print('-> ERROR:', exc)
     #       return

    # Creating a common colorbar.
    cbar_x, cbar_y, size_x, size_y = get_cbar_position(axarr) # get extent
    cbar_ax = fig.add_axes([cbar_x, cbar_y, size_x, size_y])
    scalarMap = cmx.ScalarMappable(norm=cmap_norm, cmap=cmap_key)
    cbar = fig.colorbar(scalarMap, cax=cbar_ax)
    cbar.set_label(units)
    # Creating the title from the filename
    fig.suptitle(title, fontsize='x-large')
    #fig.tight_layout()
    # Save the title
    fig.savefig(output_filename, dpi=150)
    plt.close()


def plot_grid(dataArray, output_filename, title, units, plot_type='contourf',):