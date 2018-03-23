class Matrix
  ## Wrapper class Matrix makes it easy to work with 2D lists

  attr_reader :rows, :cols

  #Fill data (leave col as 0 unless planning on setting vals)
  def initialize(rows, cols)
    @rows = rows
    @cols = cols
    @data = []
    for i in (0...rows)
      @data[i] = []
      for j in (0...cols)
        @data[i][j] = 0
      end
    end
  end

  def reset(rows, cols)
    @rows = rows
    @cols = cols
    @data = []
    for i in (0...rows)
      @data[i] = []
      for j in (0...cols)
        @data[i][j] = 0
      end
    end
  end

  # For compatablity
  def to_str to_s; end

  #Neatly represent data
  def to_s
    ret = "\n" + @rows.to_s + 'x' + @cols.to_s + " matrix:\n"
    for row in @data
      for datum in row
        ret+= "nil" if !datum
        ret+= datum.to_i.to_s + "\t"
      end
      ret+= "\n"
    end
    return ret + "\n"
  end

  #Set a coord
  def set(row, col, val)
    return nil if row >= @rows || col >= @cols || row < 0 || col < 0
    @data[row][col] = val
  end

  # Get a number
  def get(row, col)
    return nil if row >= @rows || row < 0 || col >= @cols || col < 0
    return @data[row][col]
  end

  def get_row(row)
    return nil if row >= @rows || row < 0
    return @data[row]
  end

  def get_col(col)
    return nil if col >= @cols || col < 0
    ret = []
    for i in (0...@rows)
      ret.push(@data[i][col])
    end
    return ret
  end

  #Add a collumn. Data is a list of what should be entered
  def add_col(data)
    for i in (0...@rows)
      @data[i][@cols] = data[i]
    end
    @cols += 1
  end

end
