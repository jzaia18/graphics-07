require './MatrixUtils.rb'

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
  def self.box(x, y, z, width, height, depth)
    add_edge(x, y, z, x + width, y, z)
    add_edge(x, y - height, z, x + width, y - height, z)
    add_edge(x + width, y, z, x + width, y - height, z)
    add_edge(x, y, z, x, y - height, z)
    add_edge(x, y, z - depth, x + width, y, z - depth)
    add_edge(x, y - height, z - depth, x + width, y - height, z - depth)
    add_edge(x + width, y, z - depth, x + width, y - height, z - depth)
    add_edge(x, y, z - depth, x, y - height, z - depth)
    add_edge(x, y, z, x, y, z - depth)
    add_edge(x, y - height, z, x, y - height, z - depth)
    add_edge(x + width, y, z, x + width, y, z - depth)
    add_edge(x + width, y - height, z, x + width, y - height, z - depth)
  end

  # Connects a matrix of points in a sphere-like fashion (requires gen_sphere())
  def self.sphere(cx, cy, cz, r)
    points = gen_sphere(cx, cy, cz, r)
    for i in (0...points.cols)
      p = points.get_col(i) # one single point
      add_edge(p[0], p[1], p[2], p[0], p[1], p[2])
    end
  end

  # Returns a matrix of all points on surface of a sphere (helper for sphere())
  def self.gen_sphere(cx, cy, cz, r)
    ret = Matrix.new(3, 0)
    phi = 0
    while phi < $TAU
      theta = 0
      while theta < PI
        x = r*cos(theta) + cx
        y = r*sin(theta)*cos(phi) + cy
        z = r*sin(theta)*sin(phi) + cz
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
    for i in (0...points.cols)
      p = points.get_col(i) # one single point
      add_edge(p[0], p[1], p[2], p[0], p[1], p[2])
    end
  end

  # Returns a matrix of all points on surface of a torus (helper for torus())
  def self.gen_torus(cx, cy, cz, r1, r2)
    ret = Matrix.new(3, 0)
    phi = 0
    while phi < $TAU
      theta = 0
      while theta < $TAU
        x = (r1 * cos(theta) + r2) * cos(phi) + cx
        y = r1 * sin(theta) + cy
        z = -1 * (r1 * cos(theta) + r2) * sin(phi) + cz
        ret.add_col([x, y, z])
        theta += $dt
      end
      phi += $dt
    end
    return ret
  end

  # Helper for add_edge
  def self.add_point(x, y, z)
    $EDGE_MAT.add_col([x, y, z, 1])
  end

  # Add an edge to the global edge matrix
  def self.add_edge(x0, y0, z0, x1, y1, z1)
    add_point(x0, y0, z0)
    add_point(x1, y1, z1)
  end

  # Draw the pixels in the matrix and clean it out
  def self.push_edge_matrix(edgemat: $EDGE_MAT)
    i = 0
    while i < edgemat.cols
      coord0 = edgemat.get_col(i)
      coord1 = edgemat.get_col(i + 1)
      line(coord0[0].to_i, coord0[1].to_i, coord1[0].to_i, coord1[1].to_i)
      i+=2
    end
    #edgemat.reset(4,0)
  end

end
