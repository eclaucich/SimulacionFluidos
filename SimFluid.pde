//Simulacion de fluidos en 2 dimensiones
//Fuente:
  //https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/GDC03.pdf
  
Fluid fluid;

//Control del color del fluido
float r = 0;
float g = 255;
float b = 0;
int dirR = 1;
int dirG = -1;
int dirB = 1;
float amt = .5;

//Renderizar fluido y/o campo vectorial
boolean showVectorField = false;
boolean showFluid = true;

void settings() 
{  
  size(N*SCALE, N*SCALE);
}

void setup() 
{
  colorMode(RGB, 100);
  fluid = new Fluid(.1, 0.000001, 0);
}

void mouseDragged() 
{
  fluid.AddDye(mouseX/SCALE, mouseY/SCALE, 10000);
  float amtX = mouseX - pmouseX;
  float amtY = mouseY - pmouseY;
  fluid.AddVelocity(mouseX/SCALE, mouseY/SCALE, amtX/5, amtY/5);
}

void keyPressed() 
{
  //Agregar/quitar viscosidad
  if (keyCode == UP) 
  {
    fluid.AddViscocity(0.001);
  }
  if (keyCode == DOWN) 
  {
    fluid.AddViscocity(-0.001);
  }
  
  //Controles qué se renderiza
  if (keyCode == RIGHT) 
  {
    showVectorField = !showVectorField;
  }
  if (keyCode == LEFT) 
  {
    showFluid = !showFluid;
  }
}

void draw() 
{
  background(0);
  fluid.Step();

  //Control del color
  if (r<=1)
    dirR = 1;
  if (g<=1)
    dirG = 1;
  if (b<=1)
    dirB = 1;
  if (r>=254)
    dirR = -1;
  if (g>=254)
    dirG = -1;
  if (b>=254)
    dirB = -1;

  r += dirR*amt;
  g += dirG*amt;
  b += dirB*amt;

  if (showFluid)
    fluid.RenderDye(r, g, b);
  fluid.Fade(0.3);
  if (showVectorField)
    fluid.RenderField();

}

//Dibuja una flecha con cierta rotación y longitud
void DrawArrow(int cx, int cy, int len, float angle) 
{
  pushMatrix();
  translate(cx, cy);
  rotate(angle);
  line(0, 0, len, 0);
  line(len, 0, len - 4, -4);
  line(len, 0, len - 4, 4);
  popMatrix();
}

//Calcula el ángulo necesario para rotar las flechas en el campo vectorial
float CalcAngle(float x, float y) 
{
  float a = atan(y/x);
  //'atan' trabaja en [-PI/2, PI/2] (primer y cuarto cuadrante)
  //los dos 'ifs' son para tener en cuenta el segundo y tercer cuadrante
  if (x<0 && y>0)
    a += PI;
  else if (x<0 && y<0)
    a -= PI;
  return a;
}

//Funcion que mapea una posicion en el canvas con un indice de un vector
int IX(int x, int y) 
{
  x = constrain(x, 0, N-1);
  y = constrain(y, 0, N-1);
  return x + y * N;
}
