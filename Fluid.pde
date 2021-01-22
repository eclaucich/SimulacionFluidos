final int N = 90;       //Pixeles de altura y ancho
final int SCALE = 10;   //Definicion de la imagen (con un N fijo, mayor escala -> peor definicion)
final int iter = 10;     //Cantidad de iteraciones para la resolucion del sistema de ecuaciones

class Fluid 
{
  float dt;    //paso del tiempo
  float diff;  //nivel de difusion
  float visc;  //nivel de viscosidad

  float[] s;
  float[] dye;  //

  float[] Vx;  //Velocidades actuales en 'x'
  float[] Vy;  //Velocidades actuales en 'y'

  float[] Vx0; //Velocidades iniciales en 'x'
  float[] Vy0; //Velocidades iniciales en 'y'

  Fluid(float dt, float diffusion, float viscosity)
  {
    this.dt = dt;
    this.diff = diffusion;
    this.visc = viscosity;

    this.s = new float[N*N];
    this.dye = new float[N*N];

    this.Vx = new float[N*N];
    this.Vy = new float[N*N];

    this.Vx0 = new float[N*N];
    this.Vy0 = new float[N*N];
  }

  void Step() 
  {
    //Se aplica difusion en las velocidades
    Diffuse(1, Vx0, Vx, visc, dt);
    Diffuse(2, Vy0, Vy, visc, dt);

    //Se proyectan las velocidades para mantener divergencia=0 (conservacion de la masa)
    Project(Vx0, Vy0, Vx, Vy);

    //Se aplica advección a las velocidades
    Advect(1, Vx, Vx0, Vx0, Vy0, dt);
    Advect(2, Vy, Vy0, Vx0, Vy0, dt);

    //Se proyectan las velocidades nuevamente
    Project(Vx, Vy, Vx0, Vy0);

    //Se aplica difusion y adveccion a la densidad (no requiere proyeccion)
    Diffuse(0, s, dye, diff, dt);
    Advect(0, dye, s, Vx, Vy, dt);
  }

  //Añadir viscosidad al fluido
  void AddViscocity(float amt) 
  {
    this.visc += amt;
    this.visc = constrain(this.visc, 0, 1);
  }

  //Agregar 'tinta' en cierta posicion del canvas
  //la 'tinta' hace referencia a la concentracion de fluido en una posicion
  void AddDye(int x, int y, float amt)
  {
    this.dye[IX(x, y)] += amt;
  }

  //Añadir velocidad al fluido en cierta posicion
  void AddVelocity(int x, int y, float amountX, float amountY)
  {
    int index = IX(x, y);

    this.Vx[index] += amountX;
    this.Vy[index] += amountY;
  }

  //Aplicar difusion a partir de un valor inicial (x0) y el valor actual (x)
  void Diffuse (int b, float[] x, float[] x0, float diff, float dt)
  {
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
  }

  //Proyectar velocidades, forzar a la velocidad a que conserve la masa
  //Por descomposicion de Hodge (un campo de velocidad es la suma de un campo gradiente y un campo de conservacion de masa)
  void Project(float[] velocX, float[] velocY, float[] p, float[] div)
  {
    for (int j = 1; j < N - 1; j++) {
      for (int i = 1; i < N - 1; i++) {
        div[IX(i, j)] = -0.5f*(
          velocX[IX(i+1, j)]
          -velocX[IX(i-1, j)]
          +velocY[IX(i, j+1)]
          -velocY[IX(i, j-1)]
          )/N;
        p[IX(i, j)] = 0;
      }
    }
    set_bnd(0, div); 
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);

    for (int j = 1; j < N - 1; j++) {
      for (int i = 1; i < N - 1; i++) {
        velocX[IX(i, j)] -= 0.5f * (  p[IX(i+1, j)]
          -p[IX(i-1, j)]) * N;
        velocY[IX(i, j)] -= 0.5f * (  p[IX(i, j+1)]
          -p[IX(i, j-1)]) * N;
      }
    }
    set_bnd(1, velocX);
    set_bnd(2, velocY);
  }

  //Aplica adveccion (resolvedor de la densidad)
  //Fuerza a la densidad a seguir cierto campo de velocidad
  //Se buscan hacia atras, las particulas que terminaron en el centro de cada celda y se interpola 
  void Advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY, float dt)
  {
    float i0, i1, j0, j1;

    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);

    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;

    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;

    for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
      for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
        tmp1 = dtx * velocX[IX(i, j)];
        tmp2 = dty * velocY[IX(i, j)];

        x = ifloat - tmp1; 
        y = jfloat - tmp2;

        if (x < 0.5f) x = 0.5f; 
        if (x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
        i0 = floor(x); 
        i1 = i0 + 1.0f;
        if (y < 0.5f) y = 0.5f; 
        if (y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
        j0 = floor(y);
        j1 = j0 + 1.0f; 

        s1 = x - i0; 
        s0 = 1.0f - s1; 
        t1 = y - j0; 
        t0 = 1.0f - t1;

        int i0i = int(i0);
        int i1i = int(i1);
        int j0i = int(j0);
        int j1i = int(j1);

        d[IX(i, j)] = 
          s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)])
          + s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
      }
    }

    set_bnd(b, d);
  }

  //Controla el comportamiento del fluido en los extremos del canvas
  //Contiene al fluido dentro del 'recipiente'
  //'b' indica la 'pared' del recipiente
  void set_bnd(int b, float[] x)
  {
    for (int i = 1; i < N - 1; i++) {
      x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
      x[IX(i, N-1)] = b == 2 ? -x[IX(i, N-2)] : x[IX(i, N-2)];
    }
    for (int j = 1; j < N - 1; j++) {
      x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
      x[IX(N-1, j)] = b == 1 ? -x[IX(N-2, j)] : x[IX(N-2, j)];
    }

    x[IX(0, 0)]     = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N-1)]   = 0.5f * (x[IX(1, N-1)] + x[IX(0, N-2)]);
    x[IX(N-1, 0)]   = 0.5f * (x[IX(N-2, 0)] + x[IX(N-1, 1)]);
    x[IX(N-1, N-1)] = 0.5f * (x[IX(N-2, N-1)]+ x[IX(N-1, N-2)]);
  }

  //Resolvedor del sistema de ecuaciones (Gauss-Seidel)
  //A mayor 'iter' mayor precision
  //No se invierte la matriz al ser muy dispersa, pocos elementos distintos de 0
  void lin_solve(int b, float[] x, float[] x0, float a, float c)
  {
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
      for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
          x[IX(i, j)] =
            (x0[IX(i, j)]
            + a*(x[IX(i+1, j)]
            +x[IX(i-1, j)]
            +x[IX(i, j+1)]
            +x[IX(i, j-1)]
            )) * cRecip;
        }
      }
      set_bnd(b, x);
    }
  }

  //Renderizar el fluido
  void RenderDye(float r, float g, float b) 
  {
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
        float x = i*SCALE;
        float y = j*SCALE;
        float d = this.dye[IX(i, j)];

        fill(r, g, b, d);
        noStroke();
        square(x, y, SCALE);
      }
    }
  }

  //Controla qué tanto se difumina el fluido
  void Fade(float amt) 
  {
    for (int i=0; i<this.dye.length; i++) {
      float d = this.dye[i];
      this.dye[i] = constrain(d-amt, 0, 255);
    }
  }

  //Renderizar el campo vectorial
  void RenderField() 
  {
    strokeWeight(1);
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
        DrawArrow(i*SCALE, j*SCALE, 8, CalcAngle(this.Vx[IX(i, j)], this.Vy[IX(i, j)]));
        stroke(255, 255, 255);
        if (i==5 && j==5)
          stroke(255, 0, 0);
      }
    }
  }
}
