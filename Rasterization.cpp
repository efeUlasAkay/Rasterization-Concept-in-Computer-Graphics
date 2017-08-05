
//Efe Ulas Akay Seyitoglu & Huseyin Beyan email: euase@kth.se & huseyinb@kth.se
#include <iostream>
#include "glm/glm.hpp"
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
using namespace std;
using glm::vec3;
using glm::vec2;
using glm::ivec2;
using glm::mat3;
struct Pixel

{
    
    int x;
    
    int y;
    
    float zinv;
    
    vec3 pos3d;
    
    //vec3 illumination;
    
};



float returnMax(float x,float y) {
    
    if(x > y)
        
        return x;
    
    return y;
    
    
    
}

vec3 currentNormal;

vec3 currentReflectance;

vec3 lightPos(0,-0.5,-0.7);

vec3 lightPower = 14.f*vec3( 1, 1, 1 );

vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );

struct Vertex

{
    
    vec3 position;
    
    //vec3 normal;
    
    //vec3 reflectance;
    
};



// ----------------------------------------------------------------------------



// GLOBAL VARIABLES







const int SCREEN_WIDTH = 500;



const int SCREEN_HEIGHT = 500;



SDL_Surface* screen;



int t;



vector<Triangle> triangles;



vec3 cameraPos( 0, 0, -3.001 );



mat3 R;



mat3 Rx;



float yaw = 0; // Yaw angle controlling camera rotation around y-axis



float xofYaw = 0;



vec2 a;



vector<Pixel> leftPixels;



vec3 myColor;



vector<Pixel> rightPixels;

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];





// ----------------------------------------------------------------------------



// FUNCTIONS





void Update();



void Draw();



void VertexShader( const Vertex& v, Pixel& p );



void Interpolate( Pixel a, Pixel b, vector<Pixel>& result );



void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color );



void PixelShader( const Pixel& p );





void ComputePolygonRows(
                        
                        const vector<Pixel>& vertexPixels,
                        
                        vector<Pixel>& leftPixels,
                        
                        vector<Pixel>& rightPixels
                        
                        );



void DrawRows(
              
              const vector<Pixel>& leftPixels,
              
              const vector<Pixel>& rightPixels
              
              );







void DrawRows( const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels ) {
    
    
    
    for(int i = 0; i < leftPixels.size(); i++) {
        
        
        
        Pixel delta;
        
        
        
        delta.x = glm::abs(leftPixels[i].x-rightPixels[i].x);
        
        delta.y = glm::abs(leftPixels[i].y-rightPixels[i].y);
        
        
        
        int edgePixels = glm::max(delta.x,delta.y) +1;
        
        
        
        vector<Pixel> edge(edgePixels);
        
        
        
        Interpolate(leftPixels[i],rightPixels[i],edge);
        
        
        
        for(int t = 0; t < edge.size(); t++) {
            
            
            
            //if(edge[t].zinv > depthBuffer[edge[t].x][edge[t].y]) {
            
            /*  depthBuffer[edge[t].x][edge[t].y] = edge[t].zinv;
             
             PutPixelSDL(screen,edge[t].x, edge[t].y, myColor);*/
            
            if(edge[t].x > 0 && edge[t].y > 0 && edge[t].y < SCREEN_HEIGHT && edge[t].x < SCREEN_WIDTH)
                
                PixelShader(edge[t]);
            
            //   }
            
            
            
            
            
        }
        
        
        
        
        
        
        
        
        
        
        
    }
    
    
    
}







void DrawPolygon( const vector<Vertex>& vertices )



{
    
    
    
    int V = vertices.size();
    
    
    
    vector<Pixel> vertexPixels( V );
    
    
    
    for( int i=0; i<V; i++ )
        
        
        
        VertexShader( vertices[i], vertexPixels[i] );
    
    
    
    vector<Pixel> leftPixels;
    
    
    
    vector<Pixel> rightPixels;
    
    
    
    ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
    
    
    
    DrawRows( leftPixels, rightPixels );
    
    
    
}



void ComputePolygonRows(
                        
                        const vector<Pixel>& vertexPixels,
                        
                        vector<Pixel>& leftPixels,
                        
                        vector<Pixel>& rightPixels
                        
                        )



{
    
    
    
    int max = -10000000;
    
    
    
    int min = 10000000;
    
    
    
    for(int i = 0; i < vertexPixels.size();i++) {
        
        
        
        if(vertexPixels[i].y > max) {
            
            
            
            
            
            
            
            max = vertexPixels[i].y;
            
            
            
        }
        
        
        
        
        
        
        
    }
    
    
    
    for(int i = 0; i < vertexPixels.size();i++) {
        
        
        
        if(vertexPixels[i].y < min) {
            
            
            
            min = vertexPixels[i].y;
            
            
            
        }
        
        
        
    }
    
    
    
    
    
    
    
    int numberOfRowsOccupied = max - min + 1;
    
    
    
    
    
    
    
    
    
    
    
    leftPixels.resize(numberOfRowsOccupied);
    
    
    
    
    
    
    
    rightPixels.resize(numberOfRowsOccupied);
    
    
    
    for( int i = 0; i<numberOfRowsOccupied; ++i )
        
        
        
    {
        
        
        
        leftPixels[i].x = +numeric_limits<int>::max();
        
        
        
        leftPixels[i].y = i + min;
        
        
        
        rightPixels[i].x = numeric_limits<int>::min();
        
        
        
        rightPixels[i].y = i + min;
        
        
        
    }
    
    
    
    
    
    
    
    for(int i = 0; i < vertexPixels.size(); i++) {
        
        
        
        int j = (i+1)%vertexPixels.size();
        
        
        
        Pixel delta;
        
        
        
        delta.x = glm::abs(vertexPixels[i].x-vertexPixels[j].x);
        
        delta.y = glm::abs(vertexPixels[i].y-vertexPixels[j].y);
        
        size_t edgePixels = glm::max(delta.x,delta.y) +1;
        
        
        
        vector<Pixel> edge(edgePixels);
        
        
        
        Interpolate(vertexPixels[i],vertexPixels[j],edge);
        
        
        
        for(int k=0;k<edge.size();k++){
            
            
            
            
            
            
            
            int t =edge[k].y - min;
            
            
            
            if(t < 0 || t > max)
                
                
                
                continue;
            
            
            
            if(edge[k].x < leftPixels[t].x ){
                
                
                
                leftPixels[t].x = edge[k].x;
                
                leftPixels[t].zinv = edge[k].zinv;
                
                leftPixels[t].pos3d = edge[k].pos3d;
                
            }
            
            
            
            if(edge[k].x > rightPixels[t].x ) {
                
                
                
                rightPixels[t].x = edge[k].x;
                
                rightPixels[t].zinv = edge[k].zinv;
                
                rightPixels[t].pos3d = edge[k].pos3d;
                
            }
            
        }
        
        
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    // 1. Find max and min y-value of the polygon
    
    
    
    // and compute the number of rows it occupies.
    
    
    
    // 2. Resize leftPixels and rightPixels
    
    
    
    // so that they have an element for each row.
    
    
    
    // 3. Initialize the x-coordinates in leftPixels
    
    
    
    // to some really large value and the x-coordinates
    
    
    
    // in rightPixels to some really small value.
    
    
    
    // 4. Loop through all edges of the polygon and use
    
    
    
    // linear interpolation to find the x-coordinate for
    
    
    
    // each row it occupies. Update the corresponding
    
    
    
    // values in rightPixels and leftPixels.
    
    
    
}







/*void DrawPolygonEdges( const vector<vec3>& vertices )
 
 
 
 {
 
 
 
 int V = vertices.size();
 
 
 
 // Transform each vertex from 3D world position to 2D image position:
 
 
 
 vector<ivec2> projectedVertices( V );
 
 
 
 for( int i=0; i<V; ++i )
 
 
 
 {
 
 
 
 VertexShader( vertices[i], projectedVertices[i] );
 
 
 
 }
 
 
 
 // Loop over all vertices and draw the edge from it to the next vertex:
 
 
 
 for( int i=0; i<V; ++i )
 
 
 
 {
 
 
 
 int j = (i+1)%V; // The next vertex
 
 
 
 vec3 color( 1, 1, 1 );
 
 
 
 DrawLineSDL( screen, projectedVertices[i], projectedVertices[j],
 
 
 
 color );
 
 
 
 }
 
 
 
 }
 
 */



/*
 
 void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color ) {
 
 
 
 
 
 
 
 ivec2 delta = glm::abs( a - b );
 
 
 
 int pixels = glm::max( delta.x, delta.y ) + 1;
 
 
 
 
 
 
 
 vector<ivec2> line( pixels );
 
 
 
 Interpolate( a, b, line );
 
 
 
 for(int i = 0; i < line.size(); i++) {
 
 
 
 PutPixelSDL(surface, line[i].x, line[i].y, color);
 
 
 
 }
 
 
 
 
 
 
 
 }
 
 */



void Interpolatee( vec3 a, vec3 b, vector<vec3>& result )

{
    
    int N = result.size();
    
    vec3 step = vec3(b-a) / float(max(N-1,1));
    
    vec3 current( a );
    
    for( int i=0; i<N; ++i )
        
    {
        
        result[i] = (current);
        
        current += step;
        
    }
    
}





void Interpolate( Pixel a, Pixel b, vector<Pixel>& result )



{
    
    int N = result.size();
    
    
    
    ivec2 temp1(a.x,a.y);
    
    ivec2 temp2(b.x,b.y);
    
    vec2 step = vec2(temp2-temp1) / float(max(N-1,1));
    
    float z = float(b.zinv-a.zinv)/float(max(N-1,1));
    
    //vec3 illstep = (b.illumination-a.illumination)/float(max(N-1,1));
    
    //vec3 curIlm = a.illumination;
    
    //vec3 posStep = (b.pos3d-a.pos3d)/float(max(N-1,1));
    
    //vec3 curPos = a.pos3d;
    
    
    
    vector<vec3> positions(result.size());
    
    Interpolatee(a.pos3d*a.zinv, b.pos3d*b.zinv, positions);
    
    
    
    
    
    vec2 current( temp1 );
    
    for( int i=0; i<N; ++i )
        
        
        
    {
        
        result[i].x = current.x;
        
        result[i].y = current.y;
        
        result[i].zinv = a.zinv + (i+1)*z;
        
        //result[i].illumination = curIlm;
        
        //curIlm += illstep;
        
        result[i].pos3d = positions[i]/result[i].zinv;
        
        //curPos += posStep;
        
        current += step;
        
        
        
    }
    
    
    
}











void VertexShader( const Vertex& v, Pixel& p ) {
    
    
    
    float focalPoint =  SCREEN_HEIGHT;
    
    
    
    vec3 temp = (v.position - cameraPos)*R*Rx;
    
    
    
    p.x = (float)((float)(focalPoint*temp.x)/temp.z) + SCREEN_WIDTH/2;
    
    
    
    p.y = (float)((float)(focalPoint*temp.y)/temp.z) + SCREEN_HEIGHT/2;
    
    
    
    p.zinv = ((float)(1.0/temp.z));
    
    p.pos3d = v.position;
    
    /*   vec3 normal = v.normal;
     
     float radius = glm::distance(lightPos,v.position);
     
     vec3 r = glm::normalize(lightPos-v.position);
     
     
     
     
     
     float surfaceArea = 4*3.1415f*radius*radius;
     
     vec3 d = (lightPower/surfaceArea)*returnMax(glm::dot(r,normal),0);
     
     p.illumination = (d + indirectLightPowerPerArea)*v.reflectance; */
    
}







int main( int argc, char* argv[] )



{
    
    LoadTestModel( triangles );
    
    
    
    screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
    
    
    
    t = SDL_GetTicks();
    
    
    
    while( NoQuitMessageSDL() )
        
        
        
    {
        
        
        
        Update();
        
        
        
        Draw();
        
        
        
    }
    
    
    
    
    
    
    
    SDL_SaveBMP( screen, "screenshot.bmp" );
    
    
    
    return 0;
    
    
    
}







void Update()



{
    
    
    
    // Compute frame time:
    
    
    
    int t2 = SDL_GetTicks();
    
    
    
    float dt = float(t2-t);
    
    
    
    t = t2;
    
    
    
    cout << "Render time: " << dt << " ms." << endl;
    
    
    
    Uint8* keystate = SDL_GetKeyState(0);
    
    
    
    if( keystate[SDLK_UP] )
        
        
        
        cameraPos.z +=0.01;
    
    
    
    
    
    
    
    if( keystate[SDLK_DOWN] )
        
        
        
        cameraPos.z -=0.01;
    
    
    
    
    
    
    
    if( keystate[SDLK_LEFT] )
        
        
        
        
        
        
        
    {
        
        
        
        
        
        
        
        // Move camera to the left
        
        
        
        
        
        
        
        yaw = yaw - 0.01;
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    }
    
    
    
    
    
    
    
    if( keystate[SDLK_RIGHT] )
        
        
        
        
        
        
        
    {
        
        
        
        
        
        
        
        // Move camera to the right
        
        
        
        
        
        
        
        yaw = yaw + 0.01;
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    if( keystate[SDLK_RSHIFT] )
        
        
        
        xofYaw -= 0.01;;
    
    
    
    
    
    
    
    if( keystate[SDLK_RCTRL] )
        
        
        
        xofYaw += 0.01;
    
    
    
    if( keystate[SDLK_w] )
        
        lightPos.z += 0.1;
    
    if( keystate[SDLK_s] )
        
        lightPos.z -= 0.1;
    
    if( keystate[SDLK_d] )
        
        lightPos.x -= 0.1;
    
    if( keystate[SDLK_a] )
        
        lightPos.x += 0.1;
    
    if( keystate[SDLK_q] )
        
        lightPos.y -= 0.1;
    
    if( keystate[SDLK_e] )
        
        lightPos.y += 0.1;
    
    
    
    
    
    Rx = mat3(1, 0, 0, 0,cos(xofYaw), -sin(xofYaw),0,sin(xofYaw), cos(xofYaw));
    
    
    
    R = mat3(cos(yaw), 0, sin(yaw), 0, 1, 0,-sin(yaw), 0, cos(yaw));
    
    
    
    
    
    
    
}







void Draw()



{
    
    
    
    for( int y=0; y<SCREEN_HEIGHT; ++y )
        
        for( int x=0; x<SCREEN_WIDTH; ++x )
            
            depthBuffer[y][x] = 0;
    
    
    
    
    
    SDL_FillRect( screen, 0, 0 );
    
    
    
    if( SDL_MUSTLOCK(screen) )
        
        
        
        SDL_LockSurface(screen);
    
    
    
    for( int i=0; i<triangles.size(); ++i )
        
        
        
    {
        
        
        
        myColor = triangles[i].color;
        
        
        
        vector<Vertex> vertices(3);
        
        
        
        vertices[0].position = triangles[i].v0;
        
        
        
        vertices[1].position = triangles[i].v1;
        
        
        
        vertices[2].position = triangles[i].v2;
        
        
        
        currentNormal = triangles[i].normal;
        
        currentReflectance = triangles[i].color;
        
        
        
        /*vertices[0].normal = triangles[i].normal;
         
         
         
         vertices[1].normal = triangles[i].normal;
         
         
         
         vertices[2].normal = triangles[i].normal;
         
         
         
         vertices[0].reflectance = triangles[i].color;
         
         
         
         vertices[1].reflectance = triangles[i].color;
         
         
         
         vertices[2].reflectance = triangles[i].color;*/
        
        
        
        DrawPolygon( vertices );
        
        
        
    }
    
    
    
    
    
    
    
    if ( SDL_MUSTLOCK(screen) )
        
        
        
        SDL_UnlockSurface(screen);
    
    
    
    SDL_UpdateRect( screen, 0, 0, 0, 0 );
    
    
    
    
    
    
    
}



void PixelShader( const Pixel& p )

{
    
    int x = p.x;
    
    int y = p.y;
    
    
    
    if( p.zinv > depthBuffer[y][x] ) {
        
        depthBuffer[y][x] = p.zinv;
        
        
        
        vec3 normal = currentNormal;
        
        float radius = glm::distance(lightPos,p.pos3d);
        
        vec3 r = glm::normalize(lightPos-p.pos3d);
        
        
        
        
        
        float surfaceArea = 4*3.1415f*radius*radius;
        
        vec3 d = (lightPower/surfaceArea)*returnMax(glm::dot(r,normal),0);
        
        vec3 illumination = (d + indirectLightPowerPerArea)*currentReflectance;
        
        
        
        PutPixelSDL( screen, x, y, illumination );
        
    }
    
    
    
    
    
}
