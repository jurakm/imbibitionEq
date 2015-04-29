// Mesh generation function -- test program
//
#include <iostream>
#include <sstream>
#include "mgf.hh"

#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>


// kompilacija
//  g++ -std=c++11 -g mgf.cc  -lsfml-graphics -lsfml-window -lsfml-system
int main(int argc, char** argv){

 
 
  double L=5, eps = 0.1, q = 0.2, sigma = 1, delta = 1;
  if(argc > 1){ eps = std::atof(argv[1]); if(eps== 0.0) throw; }
  if(argc > 2){ q   = std::atof(argv[2]); if(q== 0.0) throw; }
  if(argc > 3){ sigma = std::atof(argv[3]); if(sigma== 0.0) throw; }
  if(argc > 4){ delta = std::atof(argv[4]); if(delta== 0.0) throw; }

  
  MGF  d1(L, eps, q, sigma, delta);
  std::cout << " tau = " << d1.tau() << ", no iter = " << d1.iter() << std::endl;

  std::vector<double> x2;
  int N = 50;
//  d1.one_side_interval(L, N, x);
  d1.double_side_interval(N, x2);
  
//  std::cout << "One side interval\n";
//  for(auto & p : x) std::cout << p << "\n";
//
//  std::cout << "Double side interval\n";
  int j=0;
 // for(auto & p : x2){  std::cout << j << ". " << p << "\n"; j++;}

//  std::ofstream  out("MGF-out.txt");
  j=0;
  N=x2.size();
  std::vector<sf::Vertex> point(N);

  // sequence of points for drawing 
  for(auto & p : x2){ 
    point[j].position = sf::Vector2f(static_cast<double>(j)/(N-1), p);
    point[j].color = sf::Color::Red;
 //   out << static_cast<double>(j)/(N-1) <<" " << p << "\n"; 
    j++; 
  }
//  out.close();
  // Find min and max of each coordinate
  auto minmax1 = std::minmax_element(point.begin(),point.end(), 
      [](sf::Vertex const & lhs, sf::Vertex const & rhs){return lhs.position.x < rhs.position.x; });
  auto minmax2 = std::minmax_element(point.begin(),point.end(), 
      [](sf::Vertex const & lhs, sf::Vertex const & rhs){return lhs.position.y < rhs.position.y; });
 
  float  xmin = minmax1.first->position.x;
  float  xmax = minmax1.second->position.x;
  float  ymin = minmax2.first->position.y;
  float  ymax = minmax2.second->position.y;
  
  std::cout << "(" << xmin <<"," << xmax << ")-->("<< ymin << "," << ymax << ")\n";
  // SFML stuff
  sf::RenderWindow mainWindow(sf::VideoMode(800,600), "Mesh generation function");
  // Transformation to the window coordonate system
  float  X0 = 200, Y0 = 500; // Point to which (xmin, ymin) is mapped
  float  dx = xmax - xmin;
  float  dy = ymax - ymin;
  float  dX = 400, dY = 300; // extend of the data in window coordinates

  sf::Transform transform;
  // First translate, then scale 
  transform.translate(X0-(dX/dx)*xmin, Y0+(dY/dy)*ymin);
  transform.scale(dX/dx, -dY/dy);
  
  sf::RenderStates state(transform);
  // x-axis 
   float xaxismin = xmin -0.1f*dx, xaxismax = xmax +0.1f*dx;
   std::vector<sf::Vertex> xaxis(2);
   xaxis[0].position = sf::Vector2f(xaxismin, ymin);
   xaxis[0].color = sf::Color::Black;
   xaxis[1].position = sf::Vector2f(xaxismax, ymin);
   xaxis[1].color = sf::Color::Black;
  // x-axis arrow 
   std::vector<sf::Vertex> xarrow(3);
   xarrow[0].position = sf::Vector2f(xaxismax+dx/50, ymin);
   xarrow[1].position = sf::Vector2f(xaxismax, ymin+dy/70);
   xarrow[2].position = sf::Vector2f(xaxismax, ymin-dy/70);
   xarrow[0].color    = sf::Color::Black;
   xarrow[1].color    = sf::Color::Black;
   xarrow[2].color    = sf::Color::Black;
 // y-axis 
   float yaxismin = ymin -0.1f*dy, yaxismax = ymax +0.1f*dy;
   std::vector<sf::Vertex> yaxis(2);
   yaxis[0].position = sf::Vector2f(xmin, yaxismin);
   yaxis[0].color = sf::Color::Black;
   yaxis[1].position = sf::Vector2f(xmin, yaxismax);
   yaxis[1].color = sf::Color::Black;
  // y-axis arrow 
   std::vector<sf::Vertex> yarrow(3);
   yarrow[0].position = sf::Vector2f(xmin, yaxismax+dy/50);
   yarrow[1].position = sf::Vector2f(xmin+dx/70, yaxismax);
   yarrow[2].position = sf::Vector2f(xmin-dx/70, yaxismax);
   yarrow[0].color    = sf::Color::Black;
   yarrow[1].color    = sf::Color::Black;
   yarrow[2].color    = sf::Color::Black;

   // load some fonts
   sf::Font font;
   if( !font.loadFromFile("/usr/share/fonts/truetype/lyx/cmr10.ttf") ){
     std::cerr << "Cannot open fonts /usr/share/fonts/truetype/lyx/cmr10.ttf" << std::endl;
     std::exit(1);
   }
   // Coordinate labels -- the position must be transformed directly
   sf::Text xlabel("x",font,20), ylabel("y",font,20);
//   xlabel.setFont(font);
//   ylabel.setFont(font);
//   xlabel.setString("x");
//   ylabel.setString("y");
//   xlabel.setCharacterSize(12);
//   ylabel.setCharacterSize(12);
   sf::Vector2f xlabelPosition =  transform.transformPoint(xaxismax, ymin - 0.05*dy);
   sf::Vector2f ylabelPosition =  transform.transformPoint(xmin - 0.05*dx, yaxismax);
   xlabel.setPosition(xlabelPosition);
   ylabel.setPosition(ylabelPosition);
   xlabel.setColor(sf::Color::Black);
   ylabel.setColor(sf::Color::Black);

   // ticks -- set x and y tics 
   int Nxtics = 5, Nytics = 5;  // nubmer of tics
   std::vector<sf::Vertex> xtics(2*Nxtics);
   std::vector<sf::Vertex> ytics(2*Nytics);
   for(int i=0; i < xtics.size(); i+=2){
     xtics[i].position = sf::Vector2f(xmin + (i/2 +1)*dx/(Nxtics), ymin+dy/100);
     xtics[i].color = sf::Color::Black;
     xtics[i+1].position = sf::Vector2f(xmin + (i/2 +1)*dx/(Nxtics), ymin-dy/100);
     xtics[i+1].color = sf::Color::Black;
   }
   
   for(int i=0; i < ytics.size(); i+=2){
     ytics[i].position = sf::Vector2f( xmin+dx/100, ymin + (i/2 +1)*dy/Nytics);
     ytics[i].color = sf::Color::Black;
     ytics[i+1].position = sf::Vector2f(xmin-dx/100, ymin + (i/2 +1)*dy/Nytics);
     ytics[i+1].color = sf::Color::Black;
   }

   // Set numbers to tics
   std::vector<sf::Text> xlabels(Nxtics);
   
   for(int i=0; i < xlabels.size(); ++i){
     float xticvalue = xmin + (i+1)*dx/Nxtics;
     std::ostringstream buff;
     buff.setf( std::ios::fixed, std:: ios::floatfield );
     buff.precision(1);
     buff << xticvalue;
     xlabels[i].setString(buff.str());
     xlabels[i].setFont(font);
     xlabels[i].setCharacterSize(12);
     xlabels[i].setColor(sf::Color::Black);
     sf::Vector2f xlabelPosition =  transform.transformPoint(xticvalue, ymin - 0.1*dy);
     xlabels[i].setPosition(xlabelPosition);
   }
  
   std::vector<sf::Text> ylabels(Nytics);
   
   for(int i=0; i < ylabels.size(); ++i){
     float yticvalue = ymin + (i+1)*dy/Nytics;
     std::ostringstream buff;
     buff.setf( std::ios::fixed, std:: ios::floatfield );
     buff.precision(1);
     buff << yticvalue;
     ylabels[i].setString(buff.str());
     ylabels[i].setFont(font);
     ylabels[i].setCharacterSize(12);
     ylabels[i].setColor(sf::Color::Black);
     sf::Vector2f ylabelPosition =  transform.transformPoint(xmin- 0.1*dx, yticvalue);
     ylabels[i].setPosition(ylabelPosition);
   }



  while( mainWindow.isOpen() ){
    sf::Event event;
    while(mainWindow.pollEvent(event)){
      if(event.type == sf::Event::Closed) mainWindow.close();

      mainWindow.clear(sf::Color::White);
      mainWindow.draw(&point[0], point.size(),sf::LinesStrip,state);
      mainWindow.draw(&xaxis[0],2,sf::Lines, state);
      mainWindow.draw(&yaxis[0],2,sf::Lines, state);
      mainWindow.draw(&xarrow[0],3,sf::Triangles, state);
      mainWindow.draw(&yarrow[0],3,sf::Triangles, state);
      
      mainWindow.draw(xlabel);
      mainWindow.draw(ylabel);

      mainWindow.draw(&xtics[0],xtics.size(),sf::Lines, state);
      mainWindow.draw(&ytics[0],ytics.size(),sf::Lines, state);

      for(int i=0; i < xlabels.size(); ++i) mainWindow.draw(xlabels[i]);
      for(int i=0; i < ylabels.size(); ++i) mainWindow.draw(ylabels[i]);



      mainWindow.display();
    }
  }

  return 0;
}
