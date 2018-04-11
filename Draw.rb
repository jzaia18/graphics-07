require './MatrixUtils.rb'
require './VectorUtils.rb'

module Draw

  # Plot a point on GRID (from top left)
  def self.plot(x, y, r: $RC, g: $GC, b: $BC)
    y = $RESOLUTION - y
    return if x < 0 || y < 0 || x >= $RESOLUTION || y >= $RESOLUTION
    $GRID[y.to_i][x.to_i] = [r.floor, g.floor, b.floor]
  end

  # Define a line by 2 points
  def self.line(x0, y0, x1, y1, r: $RC, g: $GC, b: $BC)
    # x0 is always left of x1
    return line(x1, y1, x0, y0, r: r, g: g, b: b) if x1 < x0

    #init vars
    dy = y1-y0
    dx = x1-x0
    x = x0
    y = y0

    ## Begin actual algorithm:
    if dy >= 0 #if the line is in octants I or II
      if dy < dx #octant I
        d = 2*dy - dx
        while x < x1
          plot(x, y, r: r, g: g, b: b)
          if d > 0
            y+=1
            d-=2*dx
          end
          x+=1
          d+=2*dy
        end
        plot(x, y, r: r, g: g, b: b)

      else #octant II
        d = dy - 2*dx
        while y < y1
          plot(x, y, r: r, g: g, b: b)
          if d < 0
            x+=1
            d+=2*dy
          end
          y+=1
          d-=2*dx
        end
        plot(x, y, r: r, g: g, b: b)
      end

    else #if the line is in octants VII or VIII

      if dy.abs > dx #octant VII
        d = dy + 2*dx
        while y > y1
          plot(x, y, r: r, g: g, b: b)
          if d > 0
            x+=1
            d+=2*dy
          end
          y-=1
          d+=2*dx
        end
        plot(x, y, r: r, g: g, b: b)

      else #octant VIII
        d = 2*dy + dx
        while x < x1
          plot(x, y, r: r, g: g, b: b)
          if d < 0
            y-=1
            d+=2*dx
          end
          x+=1
          d+=2*dy
        end
        plot(x, y, r: r, g: g, b: b)
      end
    end
  end

  # Circle
  def self.circle(cx, cy, cz, rad)
    t = 0
    while (t < 1)
      add_edge(cx + rad * cos($TAU * t), cy + rad * sin($TAU * t), cz, cx + rad * cos($TAU * (t + $dt)), cy + rad * sin($TAU * (t + $dt)), cz )
      t += $dt
    end
  end

  def self.cubic(ax, bx, cx, dx, ay, by, cy, dy)
    t = 0
    while t < 1
      x0 = ax*(t**3) + bx*(t**2) + cx*t + dx
      y0 = ay*(t**3) + by*(t**2) + cy*t + dy
      t+= $dt
      x1 = ax*(t**3) + bx*(t**2) + cx*t + dx
      y1 = ay*(t**3) + by*(t**2) + cy*t + dy
      add_edge(x0, y0, 0, x1, y1, 0)
    end
  end

  def self.hermite(x0, y0, x1, y1, dx0, dy0, dx1, dy1)
    xcoors = Matrix.new(4, 1)
    xcoors.set(0, 0, x0)
    xcoors.set(1, 0, x1)
    xcoors.set(2, 0, dx0)
    xcoors.set(3, 0, dx1)
    ycoors = Matrix.new(4, 1)
    ycoors.set(0, 0, y0)
    ycoors.set(1, 0, y1)
    ycoors.set(2, 0, dy0)
    ycoors.set(3, 0, dy1)
    MatrixUtils.multiply(MatrixUtils.hermite(), xcoors)
    MatrixUtils.multiply(MatrixUtils.hermite(), ycoors)
    cubic(xcoors.get(0, 0), xcoors.get(1, 0), xcoors.get(2, 0), xcoors.get(3, 0), ycoors.get(0, 0), ycoors.get(1, 0), ycoors.get(2, 0), ycoors.get(3, 0))
  end

  def self.bezier(x0, y0, x1, y1, x2, y2, x3, y3)
    xcoors = Matrix.new(4, 1)
    xcoors.set(0, 0, x0)
    xcoors.set(1, 0, x1)
    xcoors.set(2, 0, x2)
    xcoors.set(3, 0, x3)
    ycoors = Matrix.new(4, 1)
    ycoors.set(0, 0, y0)
    ycoors.set(1, 0, y1)
    ycoors.set(2, 0, y2)
    ycoors.set(3, 0, y3)
    MatrixUtils.multiply(MatrixUtils.bezier(), xcoors)
    MatrixUtils.multiply(MatrixUtils.bezier(), ycoors)
    cubic(xcoors.get(0, 0), xcoors.get(1, 0), xcoors.get(2, 0), xcoors.get(3, 0), ycoors.get(0, 0), ycoors.get(1, 0), ycoors.get(2, 0), ycoors.get(3, 0))
  end


  # x,y,z on top left front corner
  def self.box(x0, y0, z0, width, height, depth)
    x1 = x0 + width
    y1 = y0 - height
    z1 = z0 - depth

    add_polygon(x0, y0, z0,   x1, y1, z0,   x1, y0, z0) # Front
    add_polygon(x0, y0, z0,   x0, y1, z0,   x1, y1, z0)
    add_polygon(x1, y0, z1,   x0, y1, z1,   x0, y0, z1) # Back
    add_polygon(x1, y0, z1,   x1, y1, z1,   x0, y1, z1)
    add_polygon(x1, y0, z0,   x1, y1, z1,   x1, y0, z1) # Right
    add_polygon(x1, y0, z0,   x1, y1, z0,   x0, y0, z0)
    add_polygon(x0, y0, z1,   x0, y1, z0,   x0, y0, z0) # Left
    add_polygon(x0, y0, z1,   x0, y1, z1,   x0, y1, z0)
    add_polygon(x0, y0, z1,   x1, y0, z0,   x1, y0, z1) # Top
    add_polygon(x0, y0, z1,   x0, y0, z0,   x1, y0, z0)
    add_polygon(x0, y1, z0,   x1, y1, z1,   x1, y1, z0) # Bottom
    add_polygon(x0, y1, z0,   x0, y1, z1,   x1, y1, z1)
  end

  # Connects a matrix of points in a sphere-like fashion (requires gen_sphere())
  def self.sphere(cx, cy, cz, r)
    points = gen_sphere(cx, cy, cz, r)
    i = 0
    layer_increment = (1/$dt + 1).to_i
    for i in 0...points.cols
      next if i%layer_increment == layer_increment-1
      p0 = points.get_col(i)
      p1 = points.get_col(i + 1)
      p2 = points.get_col((i + layer_increment + 1)%points.cols)
      p3 = points.get_col((i + layer_increment)%points.cols)
      add_polygon(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2])
      add_polygon(p0[0], p0[1], p0[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2])
    end
  end

  # Returns a matrix of all points on surface of a sphere (helper for sphere())
  def self.gen_sphere(cx, cy, cz, r)
    ret = Matrix.new(3, 0)
    phi = 0
    while phi < 1
      theta = 0
      _phi = phi*$TAU
      while theta < 1
        _theta = theta*PI
        x = r*cos(_theta) + cx
        y = r*sin(_theta)*cos(_phi) + cy
        z = r*sin(_theta)*sin(_phi) + cz
        ret.add_col([x, y, z])
        theta += $dt
      end
      phi += $dt
    end
    return ret
  end

  # Connects a matrix of points in a torus-like fashion (requires gen_torus())
  def self.torus(cx, cy, cz, r1, r2)
    points = gen_torus(cx, cy, cz, r1, r2)
    layer_increment = (1/$dt + 1).to_i
    for i in 0...(layer_increment-1)
      for j in 0...(layer_increment-1)
        p0 = points.get_col( i*layer_increment + j)
        p1 = points.get_col( i*layer_increment + j + 1)
        p2 = points.get_col((i*layer_increment + j + layer_increment + 1)%points.cols)
        p3 = points.get_col((i*layer_increment + j + layer_increment)%points.cols)
        add_polygon(p0[0], p0[1], p0[2], p3[0], p3[1], p3[2], p2[0], p2[1], p2[2])
        add_polygon(p0[0], p0[1], p0[2], p2[0], p2[1], p2[2], p1[0], p1[1], p1[2])
      end
    end
  end

  # Returns a matrix of all points on surface of a torus (helper for torus())
  def self.gen_torus(cx, cy, cz, r1, r2)
    ret = Matrix.new(3, 0)
    phi = 0
    while phi < $TAU
      theta = 0
      while theta <= $TAU
        x = (r1 * cos(theta) + r2) * cos(phi) + cx
        y = r1 * sin(theta) + cy
        z = -1 * (r1 * cos(theta) + r2) * sin(phi) + cz
        ret.add_col([x, y, z])
        theta += $dt * $TAU
      end
      phi += $dt * $TAU
    end
    return ret
  end

  # Helper for add_edge
  def self.add_point(x, y, z, mat)
    mat.add_col([x, y, z, 1])
  end

  # Add an edge to the global edge matrix
  def self.add_edge(x0, y0, z0, x1, y1, z1, mat: $EDGE_MAT)
    add_point(x0, y0, z0, mat)
    add_point(x1, y1, z1, mat)
  end

  # Add a trangle to the global polygon matrix
  def self.add_polygon(x0, y0, z0, x1, y1, z1, x2, y2, z2, mat: $POLY_MAT)
    add_point(x0, y0, z0, mat)
    add_point(x1, y1, z1, mat)
    add_point(x2, y2, z2, mat)
  end

  # Draw the pixels in the matrix
  def self.push_edge_matrix(edgemat: $EDGE_MAT)
    i = 0
    while i < edgemat.cols
      coord0 = edgemat.get_col(i)
      coord1 = edgemat.get_col(i + 1)
      line(coord0[0].to_i, coord0[1].to_i, coord1[0].to_i, coord1[1].to_i)
      i+=2
    end
  end

  # Draw the pixels in the matrix
  def self.push_polygon_matrix(polymat: $POLY_MAT)
    i = 0
    while i < polymat.cols
      coord0 = polymat.get_col(i)
      coord1 = polymat.get_col(i + 1)
      coord2 = polymat.get_col(i + 2)
      if (calc_normal(coord0, coord1, coord2)[2] > 0)
        line(coord0[0].to_i, coord0[1].to_i, coord1[0].to_i, coord1[1].to_i)
        line(coord1[0].to_i, coord1[1].to_i, coord2[0].to_i, coord2[1].to_i)
        line(coord2[0].to_i, coord2[1].to_i, coord0[0].to_i, coord0[1].to_i)
      end
      i+=3
    end
  end

  # Given a triangle, calculate its normal
  def self.calc_normal(p0, p1, p2)
    a = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]]
    b = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]]
    return VectorUtils.cross_product(a, b)
  end

end
