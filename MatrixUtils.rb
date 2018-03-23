include Math

module MatrixUtils

  # Create an identity matrix of side-length size
  def self.identity(size)
    ret = Matrix.new(size, size)
    for i in (0...size)
      ret.set(i, i, 1)
    end
    return ret
  end

  # Multiplies 2 matricies (dot product) and modifies second input
  def self.multiply(m1, m2, modify_second: true)
    raise "Cannot multiply these matricies due to their dimensions" if m1.cols != m2.rows
    ret = Matrix.new(m1.rows, m2.cols)
    for r in (0...ret.rows)
      for c in (0...ret.cols)
        sum = 0
        for i in (0...m1.rows) # Get the dot product and sum it
          sum += m1.get_row(r)[i] * m2.get_col(c)[i] end
        ret.set(r, c, sum)
      end
    end

    if modify_second
      for r in (0...m2.rows)
        for c in (0...m2.cols)
          m2.set(r, c, ret.get(r, c))
        end
      end
    end
    return ret
  end

  # Adds 2 matricies and modifies first input
  def self.add(m1, m2, modify_second: true)
    raise "Cannot add these matricies due to their dimensions" if m1.rows != m2.rows || m1.cols != m2.cols
    ret = Matrix.new(m1.rows, m1.cols)
    for r in (0...m1.rows)
      for c in (0...m1.cols)
        ret.set(r, c, m1.get(r, c) + m2.get(r, c))
      end
    end

    if modify_second
      for r in (0...m2.rows)
        for c in (0...m2.cols)
          m1.set(r, c, ret.get(r, c))
        end
      end
    end
    return ret
  end

  # Subtracts 2 matricies and modifies first input
  def self.subtract(m1, m2, modify_second: true)
    raise "Cannot subtract these matricies due to their dimensions" if m1.rows != m2.rows || m1.cols != m2.cols
    ret = Matrix.new(m1.rows, m1.cols)
    for r in (0...m1.rows)
      for c in (0...m1.cols)
        ret.set(r, c, m1.get(r, c) - m2.get(r, c))
      end
    end

    if modify_second
      for r in (0...m2.rows)
        for c in (0...m2.cols)
          m1.set(r, c, ret.get(r, c))
        end
      end
    end
    return ret
  end

  def self.dilation(sx, sy, sz)
    ret = identity(4);
    ret.set(0, 0, sx)
    ret.set(1, 1, sy)
    ret.set(2, 2, sz)
    return ret
  end

  def self.translation(a, b, c)
    ret = identity(4);
    ret.set(0, 3, a)
    ret.set(1, 3, b)
    ret.set(2, 3, c)
    return ret
  end

  def self.rotation(axis, theta)
    theta = theta * $TAU / 360
    ret = identity(4);

    case(axis)
    when 'x'
      ret.set(1, 1, cos(theta))
      ret.set(1, 2, -1 * sin(theta))
      ret.set(2, 1, sin(theta))
      ret.set(2, 2, cos(theta))
    when 'y'
      ret.set(0, 0, cos(theta))
      ret.set(0, 2, sin(theta))
      ret.set(2, 0, -1 * sin(theta))
      ret.set(2, 2, cos(theta))
    when 'z'
      ret.set(0, 0, cos(theta))
      ret.set(0, 1, -1 * sin(theta))
      ret.set(1, 0, sin(theta))
      ret.set(1, 1, cos(theta))
    end

    return ret
  end

  def self.hermite()
    ret = Matrix.new(4, 4)
    ret.set(0, 0, 2)
    ret.set(0, 1, -2)
    ret.set(0, 2, 1)
    ret.set(0, 3, 1)
    ret.set(1, 0, -3)
    ret.set(1, 1, 3)
    ret.set(1, 2, -2)
    ret.set(1, 3, -1)
    ret.set(2, 2, 1)
    ret.set(3, 0, 1)
    return ret
  end

  def self.bezier()
    ret = Matrix.new(4, 4)
    ret.set(0, 0, -1)
    ret.set(0, 1, 3)
    ret.set(0, 2, -3)
    ret.set(0, 3, 1)
    ret.set(1, 0, 3)
    ret.set(1, 1, -6)
    ret.set(1, 2, 3)
    ret.set(2, 0, -3)
    ret.set(2, 1, 3)
    ret.set(3, 0, 1)
    return ret
  end

end
