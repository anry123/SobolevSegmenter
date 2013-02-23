#include <stdlib.h> 
#include <vector>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "SobolevSegmenterCLP.h"

#include "itkPluginUtilities.h"


// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

//functions definition ...
#define min(a,b) (a<b)?a : b
#define max(a,b) (a<b)?b : a

void MooreBoundaryTracing(const std::vector<bool> BW, const int sz_row, const int sz_col, std::vector<float>& xi, std::vector<float>& yi);
bool isEdgePoint(const std::vector<bool> BW, const int sz_row, const int sz_col, const int row, const int col);
std::vector<float> circshift(const std::vector<float> X, const int shift);
std::vector<bool>  OutherBand(const std::vector<bool> BW, const int sz_row, const int sz_col);
std::vector<float> cconv(const std::vector<float> x1, const std::vector<float> x2);
bool isInsidePolygon(const std::vector<float> xi, const std::vector<float> yi, int x, int y);





template <class T>
int DoIt( int argc, char * argv[], T )
{
PARSE_ARGS;


	// Initialize input images and constants

    //const int iter = atoi(argv[4]); 

	//read 2 input images
	typedef   float  InputPixelType;
	typedef itk::Image< InputPixelType,  3 >   InputImageType;
	typedef itk::ImageFileReader< InputImageType  >  ReaderType;
	ReaderType::Pointer reader1 = ReaderType::New();
	ReaderType::Pointer reader2 = ReaderType::New();
	reader1->SetFileName( inputVolume.c_str() );
	reader2->SetFileName( initialMask.c_str() );
	reader1 -> Update( );
	reader2 -> Update( );
	InputImageType::Pointer image1 = reader1 -> GetOutput( );
	const InputImageType::SpacingType& sp = image1->GetSpacing();
	const InputImageType::PointType& orgn = image1->GetOrigin();
	InputImageType::Pointer image2 = reader2 -> GetOutput( );
	InputImageType::SizeType imageSize = image1->GetLargestPossibleRegion().GetSize();

	typedef itk::ImageRegionConstIterator< InputImageType > ConstIteratorType;
    ConstIteratorType in1( image1, image1->GetRequestedRegion() );
	ConstIteratorType in2( image2, image2->GetRequestedRegion() );

	int sz_row = imageSize[0], sz_col = imageSize[1];
		
	std::vector<float> I(sz_row*sz_col,0);
	std::vector<bool> BW(sz_row*sz_col,0U),mask(sz_row*sz_col,0U), BW_out(sz_row*sz_col,0U);
    std::vector<float> xi,yi;
	
	int cnt;
	for (cnt=0, in1.GoToBegin(), in2.GoToBegin(); !in1.IsAtEnd(); ++in1, ++in2, ++cnt )
    {
	    I[floor(float(cnt)/sz_col)+sz_row*(cnt % sz_col)] = in1.Get();
	    BW[floor(float(cnt)/sz_col)+sz_row*(cnt % sz_col)] = bool(in2.Get());
    }


	//int sz_row = 12, sz_col = 10;
	//std::vector<int> I(sz_row*sz_col,0);
	//std::vector<bool> BW(sz_row*sz_col,0U),mask(sz_row*sz_col,0U), BW_out(sz_row*sz_col,0U);
 //   std::vector<float> xi,yi;

	//int cnt=0;
	//int col,row;


	//for(col=0; col<sz_col; col++)
	//	for(row=0; row<sz_row; row++){
	//		if(row>=2 && col>=4 && row<=6 && col<=7) I[cnt] = 200;
	//		if(row>=3 && col>=3 && row<=7 && col<=6) BW[cnt] = 1U;
	//		cnt++;
	//	}
	
/////////////////////////////////////////////////////////////
 // HERE: I, BW, mask and iter are defined already

// Constant parameters initialization
// Kernel parameter
//float lambda = 0.1;
// The minimum force limit (in case it is smaller)
float epsilon = 1e-3;
// the 'speed' of movement, larger values result in faster convergence of the contour
// float dt = 1.0;
float mean_in,mean_out, max_force,total_length;
float tmp1,tmp2,tmp3,tmp4;
std::vector<float> arclen,cumarclen,normal_x,normal_y,force_x,force_y,force_dx,force_dy,K;
std::vector<float> xi_next, xi_prev, yi_next, yi_prev; 

    MooreBoundaryTracing(BW, sz_row, sz_col, xi, yi);

  // mask = BW;
//	mask = BW;

  // for n=1:iter
  for(int n=0; n<numiter; n++){
    // creating shifted matrices for easy usage
    //xi_next = circshift(xi, -1);
    //xi_prev = circshift(xi, 1);
    //yi_next = circshift(yi, -1);
    //yi_prev = circshift(yi, 1);
	  xi_next = circshift(xi, -1);
	  xi_prev = circshift(xi, 1);
	  yi_next = circshift(yi, -1);
	  yi_prev = circshift(yi, 1);
	// BW_out = mask; %temporary BW out
    // BW_out = conv2(float(BW_out),fspecial('gaussian',[5,5],10),'same')>0.1;
	//  BW_out = OutherBand(BW,sz_row, sz_col);
	//  BW_out = mask;
	//  for(tmp1=0;tmp1<mask.size();tmp1++) BW_out[tmp1]=!mask[tmp1];
    // mean_in = sum(sum(mask.*I))/sum(sum(mask));
    // mean_out = sum(sum(BW_out.*I))/sum(sum(BW_out))*/;


	 tmp1=0,tmp2=0,tmp3=0,tmp4=0;
 /*    for(cnt=0, in1.GoToBegin(); !in1.IsAtEnd(); ++in1,++cnt)
	 if(isInsidePolygon(yi,xi,cnt/sz_col,cnt % sz_col)){
		 tmp1 += in1.Get();
		 tmp2++;
	 }
	 else{
		 tmp3 += in1.Get();
		 tmp4++;
	 }*/

	 // detect contour bounding box
	 int minrow=sz_row,mincol=sz_col,maxrow=0,maxcol=0;
	 for(cnt=0; cnt<xi.size();cnt++){
          if (xi[cnt]>maxrow) maxrow = xi[cnt];
		  if (xi[cnt]<minrow) minrow = xi[cnt];
		  if (yi[cnt]>maxcol) maxcol = yi[cnt];
		  if (yi[cnt]<mincol) mincol = yi[cnt];
	 }
	 int row,col;
	 for(cnt=0; cnt<sz_col*sz_row; cnt++){
		 row = cnt % sz_row;
		 col = cnt/sz_row;
		 if (row>=minrow & row<=maxrow & col>=mincol & col<=maxcol)
	         if(isInsidePolygon(xi,yi,row,col)){tmp1 += I[cnt]; tmp2++;}
	         else{tmp3 += I[cnt]; tmp4++;}
		 else
		 {
         tmp3 += I[cnt];
		 tmp4++;
	     }
	 }
	 mean_in = tmp1/tmp2;
	 mean_out = tmp3/tmp4;

    // getArclen implementation - updates arclen, cumarclen and normal
    // part 1 : calculation of cumarclen and total_length
    // arclen = sqrt((xi-xi_next).^2+(yi-yi_next).^2);
    // total_length = sum(arclen);
    // cumarclen = cumsum(arclen);
	//  x normal
    // normal_x = yi_prev - yi_next;
    //  y normal
    // normal_y = xi_next - xi_prev;
	// mag = sqrt(normal_x.^2+normal_y.^2);
	// normal_x = normal_x./mag;
    // normal_y = normal_y./mag;
    total_length = 0;
	arclen = xi;
	cumarclen = xi;
	normal_x = xi;
	normal_y = yi;
	for(cnt=0;cnt<xi.size();cnt++){
        arclen[cnt] = sqrt((xi[cnt]-xi_next[cnt])*(xi[cnt]-xi_next[cnt])+(yi[cnt]-yi_next[cnt])*(yi[cnt]-yi_next[cnt]));
		total_length += arclen[cnt];
		cumarclen[cnt] = total_length;
		normal_x[cnt] = yi_prev[cnt] - yi_next[cnt];
        normal_y[cnt] = xi_next[cnt] - xi_prev[cnt];
		tmp1 = sqrt(normal_x[cnt]*normal_x[cnt]+normal_y[cnt]*normal_y[cnt]);
		normal_x[cnt] /= tmp1;
		normal_y[cnt] /= tmp1;
	}


	// chanvese implementation
    // int_xi = min(rows,max(1,round(xi)));
    // int_yi = min(cols,max(1,round(yi)));
    // Ival = I(int_yi+rows*(int_xi-1));
    // Calculate the energy to minimize
    // fval = (mean_out - mean_in).*(mean_in + mean_out - 2 * Ival);
    //  multiply with the normal in order to get the desired direction
    // force_x = normal_x.*fval;
    // force_y = normal_y.*fval;
	//  conv implementation
    //     Calculate the Kernel H1 Tilda (As can be read in the Sobolev
    //     article)
    //    K = (1 + (cumarclen.^2 - cumarclen * total_length + total_length^2/6) /...
    //        ( 2 * lambda * total_length^2)) / total_length;
	force_x = xi;
	force_y = yi;
	force_dx = xi;
	force_dy = yi;
	K = xi;
	std::vector<float> Itmp = xi;
	for(cnt=0;cnt<xi.size();cnt++){
		Itmp[cnt] = I[(min(sz_row-1,max(0,floor(0.5+xi[cnt])))) + (min(sz_col-1,max(0,floor(0.5+yi[cnt]))))*sz_row];
        force_x[cnt] = normal_x[cnt]*(mean_out-mean_in)*(mean_in+mean_out-2 * Itmp[cnt]);
		force_y[cnt] = normal_y[cnt]*(mean_out-mean_in)*(mean_in+mean_out-2 * Itmp[cnt]);
		force_dx[cnt] = force_x[cnt]*arclen[cnt];
		force_dy[cnt] = force_y[cnt]*arclen[cnt];
        K[cnt] = (1 + (cumarclen[cnt]*cumarclen[cnt] - cumarclen[cnt] * total_length + total_length*total_length/6) /
                                             ( 2.0 * lambda * total_length*total_length)) / total_length;
	}
  // force_x = cconv(force_x.*arclen,K,length(force_x));
  // force_y = cconv(force_y.*arclen,K,length(force_y));
	force_x = cconv(force_dx,K);
	force_y = cconv(force_dy,K);


	//   calculate max_force
    // max_force = max(sqrt(force_x.^2 + force_y.^2));
	max_force = 0;
	for(cnt=0;cnt<xi.size();cnt++) max_force = max(max_force, sqrt(force_x[cnt]*force_x[cnt] + force_y[cnt]*force_y[cnt]));

	// Update contour points
    // xi(:) = xi(:) + force_x .* dt ./ (max_force+epsilon);
    // yi(:) = yi(:) + force_y .* dt ./ (max_force+epsilon);
	for(cnt=0;cnt<xi.size();cnt++){
       xi[cnt] += force_x[cnt] * dt / (max_force+epsilon);
	   yi[cnt] += force_y[cnt] * dt / (max_force+epsilon);
	}


    //% create an updated roipoly
    //% BW = roipoly(BW,xi,yi);
    //mask = zeros(rows,cols);
    //for r=1:rows
    //    for c=1:cols
    //        s=numel(xi)-1;
    //        flag = false;
    //        for t=1:numel(xi)
    //            if (yi(t)< c && yi(s)>=c ||...
    //                    yi(s)< c && yi(t)>=c) &&...
    //                    (xi(t)<=r || xi(s)<=r)
    //                if (xi(t)+(c-yi(t))/(yi(s)-yi(t))*(xi(s)-xi(t))<r)
    //                    flag=~flag;
    //                end
    //            end
    //            s=t;
    //        end
    //        if flag
    //            mask(r,c) = 1;
    //        end
    //    end
    //end
    
 /*    for(int col=0;col<sz_col;col++)
		 for(int row=0;row<sz_row;row++){
			 if(isInsidePolygon(xi,yi,row,col)) mask[row+col*sz_row] = 1U;
			 else mask[row+col*sz_row] = 0U;
		 }*/
  }
  mask = BW;
     for(int col=0;col<sz_col;col++)
		 for(int row=0;row<sz_row;row++){
			 if(isInsidePolygon(xi,yi,row,col)) mask[row+col*sz_row] = 1U;
			 else mask[row+col*sz_row] = 0U;
		 }

// create output image
    typedef itk::Image< unsigned int,  3 >   OutputImageType;
	OutputImageType::Pointer imageOut = OutputImageType::New();
	OutputImageType::IndexType start;
    start[0] = 0; // first index on X
    start[1] = 0; // first index on Y
	start[2] = 0; // first index on Z
	OutputImageType::SizeType size;
    size[0] = sz_row; // size along X
    size[1] = sz_col; // size along Y
	size[2] = 1; // size along Z
	imageOut->SetOrigin( orgn );
	imageOut->SetSpacing( sp );
	OutputImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
	imageOut->SetRegions( region );
    imageOut->Allocate();

    typedef itk::ImageRegionIterator< OutputImageType > IteratorType;
    IteratorType out_iter( imageOut, imageOut->GetRequestedRegion() );
		
	for ( cnt=0, out_iter.GoToBegin(); !out_iter.IsAtEnd(); ++out_iter, ++cnt )
    {
		out_iter.Set(int(mask[floor(float(cnt)/sz_row)+sz_col*(cnt % sz_row)]));
    }

	// write image to disk
	typedef itk::ImageFileWriter< OutputImageType >  WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( outputVolume.c_str() );
	writer->SetInput( imageOut );
	writer->Update();

  return EXIT_SUCCESS;
}



void MooreBoundaryTracing(const std::vector<bool> BW, const int sz_row, const int sz_col, std::vector<float>& xi, std::vector<float>& yi)
{
   float dir[2];
   float tmp;

  //R = [ cos(pi/4) , -sin(pi/4) ; sin(pi/4) , cos(pi/4)];
   const float Rinv[4] =  { 1/sqrt(2.0),1/sqrt(2.0), -1/sqrt(2.0), 1/sqrt(2.0) }; 
   const float Rpow3[4] = { -1/sqrt(2.0),-1/sqrt(2.0), 1/sqrt(2.0), -1/sqrt(2.0) };
  
  //cnt = 1;
  //while BW(cnt)==false
  //    cnt=cnt+1;
  //end
  //[xi(1,1),yi(1,1)] = ind2sub(size(BW),cnt);
   int cnt = 0;
   while(BW[cnt]==0U) cnt++;
   xi.push_back(cnt % sz_row);
   yi.push_back(cnt/sz_row);

  //if isedge(BW,xi(1)-1,yi(1)+1)
  //    dir = [-1.0;1.0];
  //elseif isedge(BW,xi(1),yi(1)+1)
  //    dir = [0;1.0];    
  //else
  //    dir = [1.0;1.0];        
  //end
  //xi(2,1) = xi(1)+dir(1);
  //yi(2,1) = yi(1)+dir(2);
  if(isEdgePoint(BW, sz_row, sz_col, xi[0]-1,yi[0]+1)){
	  dir[0] = -1.0;
	  dir[1] = 1.0;
  }
  else 
	  if(isEdgePoint(BW, sz_row, sz_col, xi[0],yi[0]+1)){
	      dir[0] = 0.0;
	      dir[1] = 1.0;
	  }
      else{
	     dir[0] = 1.0;
	     dir[1] = 1.0;
	  }
   xi.push_back(xi[0]+dir[0]);
   yi.push_back(yi[0]+dir[1]);

  //while (xi(end,1)+dir(1))~=xi(1,1) || (yi(end,1)+dir(2))~=yi(1,1)
   while((xi.back()+dir[0]) != xi.front() || (yi.back()+dir[1]) != yi.front()){
  //    dir = round(R^3*dir);
	   tmp = dir[0];
	   dir[0] = floor(0.5+Rpow3[0]*dir[0]+Rpow3[1]*dir[1]);
	   dir[1] = floor(0.5+Rpow3[2]*tmp+Rpow3[3]*dir[1]);
  //    while ~isedge(BW,xi(end,1)+dir(1),yi(end,1)+dir(2))
	   while(!isEdgePoint(BW, sz_row, sz_col, xi.back()+dir[0],yi.back()+dir[1])){
  //        dir = round(R^-1*dir);
		     tmp = dir[0];
	         dir[0] = floor(0.5+Rinv[0]*dir[0]+Rinv[1]*dir[1]);
	         dir[1] = floor(0.5+Rinv[2]*tmp+Rinv[3]*dir[1]);
  //    end
	   }
  //    xi(end+1,1) = xi(end,1)+dir(1);
  //    yi(end+1,1) = yi(end,1)+dir(2);
        xi.push_back(xi.back()+dir[0]);
        yi.push_back(yi.back()+dir[1]);    
  //end	
   }
}



bool isEdgePoint(const std::vector<bool> BW, const int sz_row, const int sz_col, const int row, const int col) 
{
	// patch = BW(a-1:a+1,b-1:b+1);
    // flag = BW(a,b)==1 && min(patch(:))==0;
  return BW[row+col*sz_row]==1U && (BW[row-1+(col-1)*sz_row]==0U || BW[row+(col-1)*sz_row]==0U || BW[row+1+(col-1)*sz_row]==0U || BW[row+1+col*sz_row]==0U ||
	                               BW[row+1+(col+1)*sz_row]==0U || BW[row+(col+1)*sz_row]==0U || BW[row-1+(col+1)*sz_row]==0U || BW[row-1+col*sz_row]==0U );
}


std::vector<float> circshift(const std::vector<float> X, const int shift)
{
  std::vector<float> tmp;
  if(shift==1){
	  tmp.push_back(X.back());
	  for(int t = 0; t<X.size()-1; t++)  tmp.push_back(X[t]);
  }
  else{
	  for(int t = 1; t<X.size(); t++)    tmp.push_back(X[t]);
	  tmp.push_back(X.front());
  }
  return tmp;
}


std::vector<bool>  OutherBand(const std::vector<bool> BW, const int sz_row, const int sz_col)
{
std::vector<bool> C(sz_row*sz_col,0U);
bool tmp;
int row,col,sub_row,sub_col;
for(col=2;col<sz_col-2;col++)
	for(row=2;row<sz_row-2;row++){
		tmp = 0U;
	    for(sub_col=-2;sub_col<2;sub_col++)
		  for(sub_row=-2;sub_row<2;sub_row++) tmp = tmp || BW[row+sub_row+(col+sub_col)*sz_row];
		if(tmp && !BW[row+col*sz_row])  C[row+col*sz_row] = 1U;
	}
return C;
}

std::vector<float> cconv(const std::vector<float> x1, const std::vector<float> x2)
{
	std::vector<float> C(x1.size(),0);
	int len = x1.size();
	float tmp;

	for(int n=0;n<len;n++)
		for(int k=0;k<len;k++)	
			if(k>n) C[n] += x1[k]*x2[(len-k+n) % len];
			else    C[n] += x1[k]*x2[(n-k) % len];
	return C;
//
//function [c]=cconv(x1,x2)
//if length(x1)>length(x2)
//    x2=[x2, zeros(1,length(x1)-length(x2))];
//else
//    x1=[x1, zeros(1,length(x2)-length(x1))];
//end
//l=length(x1);
//n=1:l;
//y(n)=0;
//for k=1:l
//p=(mod(n-k,l))+1;
//y(n)=y(n)+ x1(k)*x2(p);
//end
//c=y;
}

bool isInsidePolygon(const std::vector<float> xi, const std::vector<float> yi, int x, int y)
{
  int   i, j=xi.size()-1 ;
  bool  flag=0U;

  for (i=0; i<xi.size(); i++) {
    if ((yi[i]< y && yi[j]>=y
    ||   yi[j]< y && yi[i]>=y)
    &&  (xi[i]<=x || xi[j]<=x)) {
      flag^=(xi[i]+(y-yi[i])/(yi[j]-yi[i])*(xi[j]-xi[i])<x); }
    j=i; }

  return flag; 

}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType(inputVolume, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch( componentType )
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0) );
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<char>(0) );
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( argc, argv, static_cast<unsigned short>(0) );
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( argc, argv, static_cast<unsigned int>(0) );
        break;
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0) );
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
