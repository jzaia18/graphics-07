require './Draw.rb'
require './Utils.rb'
require './MatrixUtils.rb'
require './Matrix.rb'

include Math

##TAU!!!!
$TAU = PI*2

# Changeable
$RESOLUTION = 500 # All images are squares
$DEBUGGING = false
$BACKGROUND_COLOR = [255, 255, 224] # [r, g, b]
$DRAW_COLOR = [255, 105, 180] # [r, g, b]
$INFILE = "script"
$OUTFILE = "image.ppm"
$TEMPFILE = "temmmmp.ppm" # Used as temp storage for displaying
$dt = 0.05 # The amount that the parametric t is incremented by on each loop

# Static
$GRID = Utils.create_grid()
$TRAN_MAT = MatrixUtils.identity(4) # Transformations matrix
$EDGE_MAT = Matrix.new(4, 0) # Edge matrix
$RC = $DRAW_COLOR[0] # Red component
$GC = $DRAW_COLOR[1]
$BC = $DRAW_COLOR[2]

##=================== MAIN ==========================
### Take in script file

if ARGV[0]
  $INFILE = ARGV[0]
end

Utils.parse_file()
