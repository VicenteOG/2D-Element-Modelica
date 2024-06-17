package Components
  connector Port
    Modelica.Units.SI.Pressure p;
    flow Modelica.Units.SI.MassFlowRate m;
    stream Modelica.Units.SI.Concentration c;
    annotation(
      Icon(graphics = {Ellipse(fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}),
      Diagram(graphics = {Ellipse(fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}));
  end Port;

  model BasicElement
    Port north annotation(
      Placement(transformation(origin = {0, 60}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}})));
    Port south annotation(
      Placement(transformation(origin = {0, -62}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {0, -100}, extent = {{-20, -20}, {20, 20}})));
    Port west annotation(
      Placement(transformation(origin = {-60, 0}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}})));
    Port east annotation(
      Placement(transformation(origin = {62, 0}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {100, 0}, extent = {{-20, -20}, {20, 20}})));
    // "Made-up" variables
    Modelica.Units.SI.Mass m_cv(start = 5*MM*V);
    Modelica.Units.SI.Concentration c(start = 5);
    Modelica.Units.SI.MassFlowRate m_dot_c(start = 0.05), m_dot_ew(start=0);
    Modelica.Units.SI.Density rho_in, rho_c(start = 1), rho_out, rho_w, rho_e;
    Modelica.Units.SI.Velocity v_in, v_out, v_w, v_e;
    Real dp_cl(unit = "Pa/m"), dp_cl_w(unit = "Pa/m");
    //Real R;
    Modelica.Units.SI.VolumeFlowRate Q(start = 1);//, Q_ew;
    parameter Modelica.Units.SI.Length l = 1, w = 1, d = 1;
    parameter Real G = 0.005;
  protected
    parameter Modelica.Units.SI.Area a_N = d*w, a_S = d*w, a_W = d*l, a_E = d*l;
    parameter Modelica.Units.SI.Volume V = l*w*d;
    parameter Modelica.Units.SI.DiffusionCoefficient diff = 1e-9;
    constant Modelica.Units.SI.MolarMass MM = 1/1000;
    constant Modelica.Units.SI.Acceleration g = 9.81;
  public
  equation
  //// ASSIGNMENTS TO "MADE-UP" VARIABLES ////
  // densities
    rho_in = inStream(south.c)*MM;
    rho_c = c*MM;
    rho_out = (north.c)*MM;
    rho_w = west.c*MM;
    rho_e = east.c*MM;
  // velocities
    v_in = south.m/rho_in/a_S;
    v_w = west.m/rho_w/a_W;
    v_e = east.m/rho_e/a_E;
  // pressure loss
    Q = south.m/rho_in;
    dp_cl = (1/d)*(rho_c*((Q/a_S)^2)/2);
    dp_cl_w = (1/d)*(rho_c*((m_dot_ew/a_E/rho_c)^2)/2);
  //// ASSIGNMENTS TO "MADE-UP" VARIABLES ////
  //// EQUATIONS ////
  // Species conservation
    der(c) = (south.m/MM/a_S + north.m/MM/a_N)/l + (west.m/MM/a_W + east.m/MM/a_E)/w + G/MM/a_E/l;
  // Mass conservation
    der(m_cv) = (south.m + north.m + west.m + east.m) + G;
  // Mass flows
    north.m = -v_out*(north.c)*a_N*MM;
    east.m = -diff*(c - east.c)/(w/2)*a_E*MM;
    west.m = -diff*(west.c - c)/(w/2)*a_W*MM;
  // MOMENTUM CONSERVATION
    l*der(m_dot_c) = (rho_in*(v_in^2)*a_S - rho_out*(v_out^2)*a_N) - (north.p*a_N - south.p*a_S) - rho_out*g*l*((a_N + a_S)/2) - dp_cl*l*((a_N + a_S)/2);
    w*der(m_dot_ew) = (rho_w*(v_w^2)*a_W - rho_e*(v_e^2)*a_E) - (east.p*a_E - west.p*a_W) - dp_cl_w*w*((a_E + a_W)/2);
  // BERNOULLI
    south.p + 0.5*(v_in^2)*rho_in = north.p + 0.5*(v_out^2)*rho_out + rho_out*g*l + dp_cl*l;
    west.p + 0.5*(v_w^2)*rho_w = east.p + 0.5*(v_e^2)*rho_e + dp_cl_w*w;
  // Pressure
    (south.p - north.p) = dp_cl*l + rho_out*g*l;
    east.p = north.p;
    
    c = (south.c + north.c)/2;
    south.c = inStream(south.c);
    annotation(
      Icon(graphics = {Rectangle(fillColor = {170, 0, 127}, fillPattern = FillPattern.CrossDiag, extent = {{-100, 100}, {100, -100}})}),
      Diagram(graphics = {Line(origin = {-0.0633566, -118.405}, points = {{0, -1.5}, {-1, 34.5}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 17), Line(origin = {0.430013, 85.3565}, points = {{0, -1.5}, {-1, 32.5}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 17), Line(origin = {-71.9575, 0.40054}, points = {{-64, 0}, {-10, 0}}, thickness = 1, arrow = {Arrow.Filled, Arrow.Filled}, arrowSize = 17), Line(origin = {147.592, -0.0928294}, points = {{-64, 0}, {-10, 0}}, thickness = 1, arrow = {Arrow.Filled, Arrow.Filled}, arrowSize = 17)}));
  end BasicElement;

  model BoundarySide
    Port west annotation(
      Placement(transformation(origin = {-60, 0}, extent = {{-40, -40}, {40, 40}}), iconTransformation(origin = {-60, 0}, extent = {{-40, -40}, {40, 40}})));
  equation
    west.m = 0;
    west.c = inStream(west.c);
    annotation(
      Icon(graphics = {Rectangle(origin = {20, 0}, fillColor = {0, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-80, 100}, {80, -100}})}));
  end BoundarySide;

  model BoundaryInlet
    Port north annotation(
      Placement(transformation(origin = {0, 60}, extent = {{-40, -40}, {40, 40}}), iconTransformation(origin = {0, 60}, extent = {{-40, -40}, {40, 40}})));
    parameter Modelica.Units.SI.Concentration C = 5;
    parameter Modelica.Units.SI.MassFlowRate m_dot = 5;
  equation
    north.m = -m_dot;
    north.c = C;
    annotation(
      Icon(graphics = {Rectangle(origin = {0, -20}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 80}, {100, -80}})}));
  end BoundaryInlet;

  model BoundaryOutlet
    Port south annotation(
      Placement(transformation(origin = {0, 60}, extent = {{-40, -40}, {40, 40}}), iconTransformation(origin = {0, -60}, extent = {{-40, -40}, {40, 40}})));
    parameter Modelica.Units.SI.Pressure P = 5e5;
  equation
    south.c = inStream(south.c);
    south.p = P;
    annotation(
      Icon(graphics = {Rectangle(origin = {0, 20}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 80}, {100, -80}})}));
  end BoundaryOutlet;

  model BasicPipe
    Port north annotation(
      Placement(transformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}})));
    Port south annotation(
      Placement(transformation(origin = {0, -100}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {0, -100}, extent = {{-20, -20}, {20, 20}})));
    Modelica.Units.SI.Mass m_cv(start = 5*MM*V);
    // "Made-up" variables
    Modelica.Units.SI.MassFlowRate m_dot(start = 5);
    Modelica.Units.SI.Density rho(start = 5*MM);
    Modelica.Units.SI.Velocity v_in, v_out;
    //, v_c;
    Real dp_cl(unit = "Pa/m");
    parameter Modelica.Units.SI.Length l = 1, w = 1, d = 1;
  protected
    parameter Modelica.Units.SI.Area a_N = 1, a_S = 1, a_W = 1, a_E = 1;
    parameter Modelica.Units.SI.Volume V = l*w*d;
    constant Modelica.Units.SI.MolarMass MM = 1/1000;
    constant Modelica.Units.SI.Acceleration g = -9.81;
  public
  equation
  //// ASSIGNMENTS TO "MADE-UP" VARIABLES ////
  // density
    rho = (south.c)*MM;
  // mass flow
    m_dot = (south.m - north.m)/2;
  // velocity
    v_in = south.m/(south.c)/a_S/MM;
  // pressure loss
    dp_cl = (1/d)*(rho*(v_in^2)/2);
  //// ASSIGNMENTS TO "MADE-UP" VARIABLES ////
  //// EQUATIONS ////
  // Mass Conservation
    der(m_cv) = (south.m + north.m);
  // MOMENTUM CONSERVATION
    l*der(m_dot) = (rho*a_S*v_in^2 - rho*a_N*v_out^2) - (north.p*a_N - south.p*a_S) - rho*g*l*((a_N + a_S)/2) - dp_cl*l*((a_N + a_S)/2);
  // PRESSURE
    (south.p - north.p) = l*dp_cl*(v_in*a_N);
  // Ports connection
    0 = north.m + south.m;
    south.c = inStream(north.c);
    north.c = inStream(south.c);
    annotation(
      Diagram,
      Icon(graphics = {Rectangle(fillColor = {0, 170, 255}, fillPattern = FillPattern.VerticalCylinder, extent = {{-40, 100}, {40, -100}})}));
  end BasicPipe;

  record SimulationInputs
    extends Modelica.Icons.Record;
    parameter Modelica.Units.SI.Concentration c = 5;
    parameter Modelica.Units.SI.MassFlowRate m = 0.5;
    parameter Modelica.Units.SI.MolarMass MM = 1/1000;
    parameter Modelica.Units.SI.Volume V = 1;
    parameter Integer n;
  end SimulationInputs;

  package Tests
    extends Modelica.Icons.ExamplesPackage;

    model Test
      extends Modelica.Icons.Example;
      BoundaryInlet boundaryInlet(m_dot = simulationInputs.m, C = simulationInputs.c) annotation(
        Placement(transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}})));
      BoundaryOutlet boundaryOutlet annotation(
        Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
      BoundarySide boundarySide_r annotation(
        Placement(transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}})));
      BoundarySide boundarySide_l annotation(
        Placement(transformation(origin = {-80, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      BasicElement basicElement(m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V), c(start = simulationInputs.c), m_dot_c(start = simulationInputs.m), G = 0, rho_c(start = simulationInputs.c*simulationInputs.MM), Q(start = simulationInputs.m/(simulationInputs.c*simulationInputs.MM))) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}})));
      SimulationInputs simulationInputs(m = 0.000025) annotation(
        Placement(transformation(origin = {64, 59}, extent = {{-32, -33}, {32, 33}})));
    equation
      connect(basicElement.west, boundarySide_l.west) annotation(
        Line(points = {{-10, 0}, {-70, 0}}));
      connect(basicElement.north, boundaryOutlet.south) annotation(
        Line(points = {{0, 10}, {0, 54}}));
      connect(basicElement.east, boundarySide_r.west) annotation(
        Line(points = {{10, 0}, {70, 0}}));
      connect(boundaryInlet.north, basicElement.south) annotation(
        Line(points = {{0, -52}, {0, -10}}));
    end Test;

    model Test2
      extends Modelica.Icons.Example;
      BoundaryInlet boundaryInlet(m_dot = simulationInputs.m, C = simulationInputs.c) annotation(
        Placement(transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}})));
      BoundaryOutlet boundaryOutlet annotation(
        Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
      BoundarySide boundarySide_r annotation(
        Placement(transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}})));
      BoundarySide boundarySide_l annotation(
        Placement(transformation(origin = {-80, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      BasicElement basicElement(G = 0, c(start = simulationInputs.c), m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V/2), w = 0.5, m_dot_c(start = simulationInputs.m/2), rho_c(start = simulationInputs.c*simulationInputs.MM), Q(start = (simulationInputs.m/2)/(simulationInputs.c*simulationInputs.MM))) annotation(
        Placement(transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
      BasicElement basicElement1(G = 0, c(start = simulationInputs.c), m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V/2), w = 0.5, m_dot_c(start = simulationInputs.m/2), rho_c(start = simulationInputs.c*simulationInputs.MM), Q(start = (simulationInputs.m/2)/(simulationInputs.c*simulationInputs.MM))) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -0)));
      SimulationInputs simulationInputs(m = 0.000025) annotation(
        Placement(transformation(origin = {67, 69}, extent = {{-29, -27}, {29, 27}})));
    equation
      connect(basicElement.north, boundaryOutlet.south) annotation(
        Line(points = {{-40, 10}, {-40, 40}, {0, 40}, {0, 54}}));
      connect(basicElement1.north, boundaryOutlet.south) annotation(
        Line(points = {{0, 10}, {0, 54}}));
      connect(boundaryInlet.north, basicElement1.south) annotation(
        Line(points = {{0, -52}, {0, -10}}));
      connect(boundaryInlet.north, basicElement.south) annotation(
        Line(points = {{0, -52}, {0, -40}, {-40, -40}, {-40, -10}}));
      connect(basicElement.west, boundarySide_l.west) annotation(
        Line(points = {{-50, 0}, {-74, 0}}));
      connect(basicElement1.west, basicElement.east) annotation(
        Line(points = {{-10, 0}, {-30, 0}}));
      connect(basicElement1.east, boundarySide_r.west) annotation(
        Line(points = {{10, 0}, {74, 0}}));
    end Test2;
    
    model Test3
      extends Modelica.Icons.Example;
      BoundaryInlet boundaryInlet(m_dot = simulationInputs.m, C = simulationInputs.c) annotation(
        Placement(transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}})));
      BoundaryOutlet boundaryOutlet annotation(
        Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
      BoundarySide boundarySide_r annotation(
        Placement(transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}})));
      BoundarySide boundarySide_l annotation(
        Placement(transformation(origin = {-80, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      BasicElement basicElement(G = 0, c(start = simulationInputs.c), m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V/simulationInputs.n), w = 1/simulationInputs.n, m_dot_c(start = simulationInputs.m/simulationInputs.n), rho_c(start = simulationInputs.c*simulationInputs.MM), Q(start = (simulationInputs.m/simulationInputs.n)/(simulationInputs.c*simulationInputs.MM))) annotation(
        Placement(transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
      BasicElement basicElement1(G = 0, c(start = simulationInputs.c), m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V/simulationInputs.n), w = 1/simulationInputs.n, m_dot_c(start = simulationInputs.m/simulationInputs.n), rho_c(start = simulationInputs.c*simulationInputs.MM), Q(start = (simulationInputs.m/simulationInputs.n)/(simulationInputs.c*simulationInputs.MM))) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -0)));
      SimulationInputs simulationInputs(m = 0.000025, n = 3) annotation(
        Placement(transformation(origin = {73, 74}, extent = {{-23, -22}, {23, 22}})));
  BasicElement basicElement11(G = 0, Q(start = (simulationInputs.m/simulationInputs.n)/(simulationInputs.c*simulationInputs.MM)), c(start = simulationInputs.c), m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V/simulationInputs.n), m_dot_c(start = simulationInputs.m/simulationInputs.n), rho_c(start = simulationInputs.c*simulationInputs.MM), w = 1/simulationInputs.n) annotation(
        Placement(transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(basicElement.north, boundaryOutlet.south) annotation(
        Line(points = {{-40, 10}, {-40, 40}, {0, 40}, {0, 54}}));
      connect(basicElement1.north, boundaryOutlet.south) annotation(
        Line(points = {{0, 10}, {0, 54}}));
      connect(boundaryInlet.north, basicElement1.south) annotation(
        Line(points = {{0, -52}, {0, -10}}));
      connect(boundaryInlet.north, basicElement.south) annotation(
        Line(points = {{0, -52}, {0, -40}, {-40, -40}, {-40, -10}}));
      connect(basicElement.west, boundarySide_l.west) annotation(
        Line(points = {{-50, 0}, {-74, 0}}));
      connect(basicElement1.west, basicElement.east) annotation(
        Line(points = {{-10, 0}, {-30, 0}}));
  connect(boundaryInlet.north, basicElement11.south) annotation(
        Line(points = {{0, -56}, {0, -40}, {40, -40}, {40, -10}}));
  connect(basicElement1.east, basicElement11.west) annotation(
        Line(points = {{10, 0}, {30, 0}}));
  connect(basicElement11.north, boundaryOutlet.south) annotation(
        Line(points = {{40, 10}, {40, 40}, {0, 40}, {0, 54}}));
  connect(basicElement11.east, boundarySide_r.west) annotation(
        Line(points = {{50, 0}, {74, 0}}));
    end Test3;

    model Test4
      extends Modelica.Icons.Example;
      BoundaryInlet boundaryInlet(m_dot = simulationInputs.m, C = simulationInputs.c) annotation(
        Placement(transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}})));
      BoundaryOutlet boundaryOutlet annotation(
        Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
      BoundarySide boundarySide_r annotation(
        Placement(transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}})));
      BoundarySide boundarySide_l annotation(
        Placement(transformation(origin = {-80, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      BasicElement basicElement(m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V/2), c(start = simulationInputs.c), m_dot_c(start = simulationInputs.m/2), G = 0, w = 0.5) annotation(
        Placement(transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}})));
      BasicElement basicElement1(m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V/2), c(start = simulationInputs.c), G = 0, m_dot_c(start = simulationInputs.m/2), w = 0.5) annotation(
        Placement(transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}})));
      BoundarySide boundarySide annotation(
        Placement(transformation(origin = {-12, 0}, extent = {{-10, -10}, {10, 10}})));
      BoundarySide boundarySide1 annotation(
        Placement(transformation(origin = {12, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Components.SimulationInputs simulationInputs annotation(
        Placement(transformation(origin = {71, 77}, extent = {{-21, -23}, {21, 23}})));
    equation
      connect(basicElement.west, boundarySide_l.west) annotation(
        Line(points = {{-50, 0}, {-70, 0}}));
      connect(basicElement.north, boundaryOutlet.south) annotation(
        Line(points = {{-40, 10}, {-40, 40}, {0, 40}, {0, 54}}));
      connect(basicElement1.east, boundarySide_r.west) annotation(
        Line(points = {{20, 0}, {70, 0}}));
      connect(basicElement1.north, boundaryOutlet.south) annotation(
        Line(points = {{40, 10}, {40, 40}, {0, 40}, {0, 54}}));
      connect(boundaryInlet.north, basicElement.south) annotation(
        Line(points = {{0, -52}, {0, -40}, {-40, -40}, {-40, -10}}));
      connect(boundaryInlet.north, basicElement1.south) annotation(
        Line(points = {{0, -52}, {0, -40}, {40, -40}, {40, -10}}));
      connect(boundarySide.west, basicElement.east) annotation(
        Line(points = {{-22, 0}, {-30, 0}}));
      connect(boundarySide1.west, basicElement1.west) annotation(
        Line(points = {{22, 0}, {30, 0}}));
    end Test4;
  end Tests;
end Components;
