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
    Modelica.Units.SI.MassFlowRate m_dot_c(start = 0.05);
    Modelica.Units.SI.Density rho_in, rho_c, rho_out;
    Modelica.Units.SI.Velocity v_in, v_c, v_out;
    Real dp_cl(unit = "Pa/m");
    
    parameter Modelica.Units.SI.Length l=1, w=1, d=1;
    parameter Real G = 0.005;
  protected
    parameter Modelica.Units.SI.Area a_N = 1, a_S = 1, a_W = 1, a_E = 1;
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
    rho_out = north.c*MM;
    // mass flow
    m_dot_c = (south.m - north.m)/2;
    // velocities
    v_c = m_dot_c/c/a_S/MM;
    v_in = south.m/inStream(south.c)/a_S/MM;
    // pressure loss
    dp_cl = (1/d)*(rho_c*(v_c^2)/2);
    //// ASSIGNMENTS TO "MADE-UP" VARIABLES ////
    
    
    //// EQUATIONS ////
    // Species conservation
    der(c) = (south.m/MM/a_S + north.m/MM/a_N)/l + (west.m/MM/a_W + east.m/MM/a_E)/w + G/MM/a_E/l;
    // Mass conservation
    der(m_cv) = (south.m + north.m + west.m + east.m) + G;
    // Mass flows
    north.m = -v_out*north.c*a_N*MM;
    east.m = -diff*(c - east.c)/(w/2)*a_E*MM;
    west.m = -diff*(west.c - c)/(w/2)*a_W*MM;
    // MOMENTUM CONSERVATION
    l*der(m_dot_c) = (rho_in*a_S*v_in^2 - rho_out*a_N*v_out^2) - (north.p*a_N - south.p*a_S) - rho_c*g*l*((a_N + a_S)/2) - dp_cl*l*((a_N + a_S)/2);
    // PRESSURE
    (south.p - north.p) = (m_dot_c^2)*l*dp_cl/((a_S^2)*(v_c^2))/(rho_c^2);
    west.p = north.p;
    east.p = north.p;
    
    //c = (east.c + west.c)/2;
    c = (south.c + north.c)/2;
    south.c = inStream(south.c);
  annotation(
      Icon(graphics = {Rectangle(fillColor = {170, 0, 127}, fillPattern = FillPattern.CrossDiag, extent = {{-100, 100}, {100, -100}})}),
 Diagram(graphics = {Line(origin = {-0.0633566, -118.405}, points = {{0, -1.5}, {(0 - 15), 32.5}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 17), Line(origin = {0.430013, 85.3565}, points = {{0, -1.5}, {(0 - 15), 32.5}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 17), Line(origin = {-71.9575, 0.40054}, points = {{-64, 0}, {-10, 0}}, thickness = 1, arrow = {Arrow.Filled, Arrow.Filled}, arrowSize = 17), Line(origin = {147.592, -0.0928294}, points = {{-64, 0}, {-10, 0}}, thickness = 1, arrow = {Arrow.Filled, Arrow.Filled}, arrowSize = 17)}));
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
    Modelica.Units.SI.Velocity v_in, v_out;//, v_c;
    Real dp_cl(unit = "Pa/m");
    
    parameter Modelica.Units.SI.Length l=1, w=1, d=1;
  protected
    parameter Modelica.Units.SI.Area a_N = 1, a_S = 1, a_W = 1, a_E = 1;
    parameter Modelica.Units.SI.Volume V = l*w*d;
    constant Modelica.Units.SI.MolarMass MM = 1/1000;
    constant Modelica.Units.SI.Acceleration g = 9.81;
    public
  equation
    //// ASSIGNMENTS TO "MADE-UP" VARIABLES ////
    // density
    rho = inStream(south.c)*MM;
    // mass flow
    m_dot = (south.m - north.m)/2;
    // velocity
    v_in = south.m/inStream(south.c)/a_S/MM;
    // pressure loss
    dp_cl = (1/d)*(rho*(v_in^2)/2);
    //// ASSIGNMENTS TO "MADE-UP" VARIABLES ////
    
    
    //// EQUATIONS ////
    // Mass Conservation
    der(m_cv) = (south.m + north.m);
    // MOMENTUM CONSERVATION
    l*der(m_dot) = (rho*a_S*v_in^2 - rho*a_N*v_out^2) - (north.p*a_N - south.p*a_S) - rho*g*l*((a_N + a_S)/2) - dp_cl*l*((a_N + a_S)/2);
    // PRESSURE
    (south.p - north.p) = (m_dot^2)*l*dp_cl/((a_S^2)*(v_in^2))/(rho^2);
    
    // Ports connection
    0 = south.m + north.m;
    south.c = inStream(north.c);
    north.c = inStream(south.c);
  annotation(
      Diagram,
      Icon(graphics = {Rectangle(fillColor = {0, 170, 255}, fillPattern = FillPattern.VerticalCylinder, extent = {{-40, 100}, {40, -100}})}));
end BasicPipe;

  record SimulationInputs
    extends Modelica.Icons.Record;
    parameter Modelica.Units.SI.Concentration c = 5;
    parameter Modelica.Units.SI.MassFlowRate m = 0.05;
    parameter Modelica.Units.SI.MolarMass MM = 1/1000;
    parameter Modelica.Units.SI.Volume V = 1;
    
  end SimulationInputs;

  package Tests
  extends Modelica.Icons.ExamplesPackage;
  
  
  model Test01
    extends Modelica.Icons.Example;
    BoundaryInlet boundaryInlet annotation(
      Placement(transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}})));
    BoundaryOutlet boundaryOutlet annotation(
      Placement(transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
    BoundaryInlet boundaryInlet1(C = 10, m_dot = 8) annotation(
      Placement(transformation(origin = {40, -62}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(boundaryInlet.north, boundaryOutlet.south) annotation(
      Line(points = {{0, -52}, {0, 34}}));
    connect(boundaryInlet1.north, boundaryOutlet.south) annotation(
      Line(points = {{40, -52}, {40, 2}, {0, 2}, {0, 34}}));
  annotation(
      Diagram(graphics = {Text(origin = {0, 78}, extent = {{-92, 19}, {92, -19}}, textString = "Example to test if the \"BoundaryOutlet\" component
\"receives\" the flow and concentration given
different mass flow and concentration comming from
the \"BoundaryInlet\" components")}));
end Test01;

  model Test02
    extends Modelica.Icons.Example;
    BoundaryInlet boundaryInlet annotation(
      Placement(transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}})));
    BoundaryOutlet boundaryOutlet annotation(
      Placement(transformation(origin = {0, 22}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
    BoundaryOutlet boundaryOutlet1 annotation(
      Placement(transformation(origin = {60, 22}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
  equation
    connect(boundaryInlet.north, boundaryOutlet.south) annotation(
      Line(points = {{0, -52}, {0, 16}}));
    connect(boundaryOutlet1.south, boundaryInlet.north) annotation(
      Line(points = {{60, 16}, {60, 0}, {0, 0}, {0, -52}}));
  annotation(
      Diagram(graphics = {Text(origin = {0, 78}, extent = {{-92, 19}, {92, -19}}, textString = "Example to test if the flow is distributed when providing
flow to two different components from one \"BoundaryInlet\"
Note: it fails! Requires a \"BasicPipe\" component that provides
a flow => \"Test03_1\" and \"Test03_2\"")}));
end Test02;
  
  model Test03_1
    extends Modelica.Icons.Example;
    BoundaryInlet boundaryInlet annotation(
      Placement(transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}})));
    BoundaryOutlet boundaryOutlet annotation(
      Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
  BasicPipe basicPipe annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}})));
  equation
    connect(boundaryInlet.north, basicPipe.south) annotation(
      Line(points = {{0, -56}, {0, -10}}));
    connect(basicPipe.north, boundaryOutlet.south) annotation(
      Line(points = {{0, 10}, {0, 54}}));
  end Test03_1;
  
  model Test03_2
    extends Modelica.Icons.Example;
    BoundaryInlet boundaryInlet annotation(
      Placement(transformation(origin = {0, -60}, extent = {{-10, -10}, {10, 10}})));
    BoundaryOutlet boundaryOutlet annotation(
      Placement(transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
    BoundaryOutlet boundaryOutlet1 annotation(
      Placement(transformation(origin = {60, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
 BasicPipe basicPipe annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}})));
 BasicPipe basicPipe1(l = 5)  annotation(
      Placement(transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}})));
  equation
 connect(boundaryInlet.north, basicPipe.south) annotation(
      Line(points = {{0, -54}, {0, -10}}));
 connect(boundaryInlet.north, basicPipe1.south) annotation(
      Line(points = {{0, -54}, {0, -40}, {60, -40}, {60, -10}}));
 connect(basicPipe.north, boundaryOutlet.south) annotation(
      Line(points = {{0, 10}, {0, 34}}));
 connect(basicPipe1.north, boundaryOutlet1.south) annotation(
      Line(points = {{60, 10}, {60, 34}}));
  annotation(
      Diagram(graphics = {Text(origin = {0, 76}, extent = {{-92, 19}, {92, -19}}, textString = "Pipes with different lengths to test if the mass flow coming
out of the inlet component is distributed in the two pipes")}));
end Test03_2;

  model Test
  extends Modelica.Icons.Example;
BoundaryInlet boundaryInlet(m_dot = simulationInputs.m, C = simulationInputs.c)  annotation(
      Placement(transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}})));
BoundaryOutlet boundaryOutlet annotation(
      Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
BoundarySide boundarySide_r annotation(
      Placement(transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}})));
BoundarySide boundarySide_l annotation(
      Placement(transformation(origin = {-80, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
BasicElement basicElement(m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V), c(start = simulationInputs.c), m_dot_c(start = simulationInputs.m))  annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}})));
 SimulationInputs simulationInputs(m = 0.05)  annotation(
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
  BoundaryInlet boundaryInlet(m_dot = 0.05)  annotation(
      Placement(transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}})));
BoundaryOutlet boundaryOutlet annotation(
      Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
  BoundarySide boundarySide_r annotation(
      Placement(transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}})));
  BoundarySide boundarySide_l annotation(
      Placement(transformation(origin = {-80, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  BasicElement basicElement(m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V), c(start = simulationInputs.c), G = 0)  annotation(
      Placement(transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
BasicElement basicElement1(m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V), c(start = simulationInputs.c), G = 0)  annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -0)));
 SimulationInputs simulationInputs(m = 0.5)  annotation(
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
  BoundaryInlet boundaryInlet(m_dot = simulationInputs.m, C = simulationInputs.c)  annotation(
      Placement(transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}})));
  BoundaryOutlet boundaryOutlet annotation(
      Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
  BoundarySide boundarySide_r annotation(
      Placement(transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}})));
  BoundarySide boundarySide_l annotation(
      Placement(transformation(origin = {-80, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  BasicElement basicElement(m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V, fixed = true), c(start = simulationInputs.c, fixed = true))  annotation(
      Placement(transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}})));
  BasicElement basicElement1(m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V, fixed = true), l = 0.15, c(start = simulationInputs.c, fixed = true))  annotation(
      Placement(transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}})));
BoundarySide boundarySide annotation(
      Placement(transformation(origin = {-12, 0}, extent = {{-10, -10}, {10, 10}})));
BoundarySide boundarySide1 annotation(
      Placement(transformation(origin = {12, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
 Components.SimulationInputs simulationInputs(m = 1)  annotation(
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
  end Test3;
    
    model Test4
    extends Modelica.Icons.Example;
    BoundaryInlet boundaryInlet(m_dot = simulationInputs.m, C = simulationInputs.c)  annotation(
        Placement(transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}})));
    BoundaryOutlet boundaryOutlet annotation(
        Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
    BoundarySide boundarySide_r annotation(
        Placement(transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}})));
    BoundarySide boundarySide_l annotation(
        Placement(transformation(origin = {-80, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    BasicElement basicElement(m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V), c(start = simulationInputs.c))  annotation(
        Placement(transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}})));
    BasicElement basicElement1(m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V), c(start = simulationInputs.c))  annotation(
        Placement(transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}})));
    BoundarySide boundarySide annotation(
        Placement(transformation(origin = {-12, 0}, extent = {{-10, -10}, {10, 10}})));
    BoundarySide boundarySide1 annotation(
        Placement(transformation(origin = {12, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    Components.SimulationInputs simulationInputs annotation(
          Placement(transformation(origin = {71, 77}, extent = {{-21, -23}, {21, 23}})));
 BasicPipe basicPipe(m_cv(fixed = true), rho(fixed = true))  annotation(
        Placement(transformation(origin = {-40, -30}, extent = {{-10, -10}, {10, 10}})));
 BasicPipe basicPipe1(m_cv(fixed = true), rho(fixed = true))  annotation(
        Placement(transformation(origin = {40, -30}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(basicElement.west, boundarySide_l.west) annotation(
        Line(points = {{-50, 0}, {-70, 0}}));
      connect(basicElement.north, boundaryOutlet.south) annotation(
        Line(points = {{-40, 10}, {-40, 40}, {0, 40}, {0, 54}}));
      connect(basicElement1.east, boundarySide_r.west) annotation(
        Line(points = {{20, 0}, {70, 0}}));
    connect(basicElement1.north, boundaryOutlet.south) annotation(
        Line(points = {{40, 10}, {40, 40}, {0, 40}, {0, 54}}));
    connect(boundarySide.west, basicElement.east) annotation(
        Line(points = {{-22, 0}, {-30, 0}}));
    connect(boundarySide1.west, basicElement1.west) annotation(
        Line(points = {{22, 0}, {30, 0}}));
 connect(boundaryInlet.north, basicPipe.south) annotation(
        Line(points = {{0, -56}, {0, -48}, {-40, -48}, {-40, -40}}));
 connect(basicPipe.north, basicElement.south) annotation(
        Line(points = {{-40, -20}, {-40, -10}}));
 connect(boundaryInlet.north, basicPipe1.south) annotation(
        Line(points = {{0, -56}, {0, -48}, {40, -48}, {40, -40}}));
 connect(basicPipe1.north, basicElement1.south) annotation(
        Line(points = {{40, -20}, {40, -10}}));
    end Test4;
    
    model Test5
    extends Modelica.Icons.Example;
    BoundaryInlet boundaryInlet(m_dot = simulationInputs.m, C = simulationInputs.c)  annotation(
        Placement(transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}})));
    BoundaryOutlet boundaryOutlet annotation(
        Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
    BoundarySide boundarySide_r annotation(
        Placement(transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}})));
    BoundarySide boundarySide_l annotation(
        Placement(transformation(origin = {-80, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    BasicElement basicElement(m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V), c(start = simulationInputs.c))  annotation(
        Placement(transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}})));
    BasicElement basicElement1(m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V), l = 1, c(start = simulationInputs.c))  annotation(
        Placement(transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}})));
    BoundarySide boundarySide annotation(
        Placement(transformation(origin = {-12, 0}, extent = {{-10, -10}, {10, 10}})));
    BoundarySide boundarySide1 annotation(
        Placement(transformation(origin = {12, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    Components.SimulationInputs simulationInputs annotation(
          Placement(transformation(origin = {71, 77}, extent = {{-21, -23}, {21, 23}})));
    BasicPipe basicPipe(m_cv(fixed = true), rho(fixed = true))  annotation(
        Placement(transformation(origin = {-40, 30}, extent = {{-10, -10}, {10, 10}})));
    BasicPipe basicPipe1(m_cv(fixed = true), rho(fixed = true))  annotation(
        Placement(transformation(origin = {40, 30}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(basicElement.west, boundarySide_l.west) annotation(
        Line(points = {{-50, 0}, {-70, 0}}));
      connect(basicElement1.east, boundarySide_r.west) annotation(
        Line(points = {{20, 0}, {70, 0}}));
    connect(boundarySide.west, basicElement.east) annotation(
        Line(points = {{-22, 0}, {-30, 0}}));
    connect(boundarySide1.west, basicElement1.west) annotation(
        Line(points = {{22, 0}, {30, 0}}));
 connect(boundaryInlet.north, basicElement.south) annotation(
        Line(points = {{0, -56}, {0, -28}, {-40, -28}, {-40, -10}}));
 connect(boundaryInlet.north, basicElement1.south) annotation(
        Line(points = {{0, -56}, {0, -28}, {40, -28}, {40, -10}}));
 connect(basicElement.north, basicPipe.south) annotation(
        Line(points = {{-40, 10}, {-40, 20}}));
 connect(basicElement1.north, basicPipe1.south) annotation(
        Line(points = {{40, 10}, {40, 20}}));
 connect(basicPipe.north, boundaryOutlet.south) annotation(
        Line(points = {{-40, 40}, {-40, 46}, {0, 46}, {0, 54}}));
 connect(basicPipe1.north, boundaryOutlet.south) annotation(
        Line(points = {{40, 40}, {40, 46}, {0, 46}, {0, 54}}));
    end Test5;
    
    model Test6
    extends Modelica.Icons.Example;
    BoundaryInlet boundaryInlet(m_dot = simulationInputs.m, C = simulationInputs.c)  annotation(
        Placement(transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}})));
    BoundaryOutlet boundaryOutlet annotation(
        Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
    BoundarySide boundarySide_r annotation(
        Placement(transformation(origin = {80, -34}, extent = {{-10, -10}, {10, 10}})));
    BoundarySide boundarySide_l annotation(
        Placement(transformation(origin = {-80, -34}, extent = {{10, -10}, {-10, 10}})));
    BasicElement basicElement(m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V), c(start = simulationInputs.c), m_dot_c(start = simulationInputs.m), l = 0.5)  annotation(
        Placement(transformation(origin = {0, -34}, extent = {{-10, -10}, {10, 10}})));
    SimulationInputs simulationInputs(m = 0.05)  annotation(
          Placement(transformation(origin = {64, 59}, extent = {{-32, -33}, {32, 33}})));
 Components.BasicElement basicElement1(c(start = simulationInputs.c), m_cv(start = simulationInputs.c*simulationInputs.MM*simulationInputs.V), m_dot_c(start = simulationInputs.m), l = 0.5) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}})));
 Components.BoundarySide boundarySide_l1 annotation(
        Placement(transformation(origin = {-80, 0}, extent = {{10, -10}, {-10, 10}})));
 Components.BoundarySide boundarySide_r1 annotation(
        Placement(transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}})));
    equation
    connect(basicElement.west, boundarySide_l.west) annotation(
        Line(points = {{-10, -34}, {-74, -34}}));
    connect(basicElement.east, boundarySide_r.west) annotation(
        Line(points = {{10, -34}, {74, -34}}));
    connect(boundaryInlet.north, basicElement.south) annotation(
        Line(points = {{0, -52}, {0, -44}}));
 connect(boundarySide_l1.west, basicElement1.west) annotation(
        Line(points = {{-74, 0}, {-10, 0}}));
 connect(basicElement1.east, boundarySide_r1.west) annotation(
        Line(points = {{10, 0}, {74, 0}}));
 connect(basicElement.north, basicElement1.south) annotation(
        Line(points = {{0, -24}, {0, -10}}));
 connect(basicElement1.north, boundaryOutlet.south) annotation(
        Line(points = {{0, 10}, {0, 54}}));
    end Test6;
  
  end Tests;

end Components;
