//******************************************************************************
// Copyright (C) 2016-2023 University of Oklahoma Board of Trustees.
//******************************************************************************
// Last modified: Fri Mar 10 18:48:57 2023 by Chris Weaver
//******************************************************************************
// Major Modification History:
//
// 20160209 [weaver]:	Original file.
// 20190203 [weaver]:	Updated to JOGL 2.3.2 and cleaned up.
// 20190227 [weaver]:	Updated to use model and asynchronous event handling.
// 20190318 [weaver]:	Modified for homework04.
// 20210320 [weaver]:	Added basic keyboard hints to drawMode().
// 20220311 [weaver]:	Improved hint wording in updatePointWithReflection().
// 20230310 [weaver]:	Improved TODO guidance especially for members to add.
//
//******************************************************************************
// Notes:
//
//******************************************************************************

package edu.ou.cs.cg.assignment.homework04;

//import java.lang.*;
import java.awt.*;
import java.awt.geom.Point2D;
import java.text.DecimalFormat;
import java.util.*;
import com.jogamp.opengl.*;
import com.jogamp.opengl.awt.GLJPanel;
import com.jogamp.opengl.glu.*;
import com.jogamp.opengl.util.FPSAnimator;
import com.jogamp.opengl.util.awt.TextRenderer;
import edu.ou.cs.cg.utilities.Utilities;

//******************************************************************************

/**
 * The <CODE>View</CODE> class.<P>
 *
 * @author  Chris Weaver
 * @version %I%, %G%
 */
public final class View
	implements GLEventListener
{
	//**********************************************************************
	// Private Class Members
	//**********************************************************************

	private static final int			DEFAULT_FRAMES_PER_SECOND = 60;
	private static final DecimalFormat	FORMAT = new DecimalFormat("0.000");

	//**********************************************************************
	// Public Class Members
	//**********************************************************************

	public static final int			MIN_SIDES = 3;
	public static final int			MAX_SIDES = 12;

	//**********************************************************************
	// Private Members
	//**********************************************************************

	// State (internal) variables
	private final GLJPanel				canvas;
	private int						w;			// Canvas width
	private int						h;			// Canvas height

	private TextRenderer				renderer;

	private final FPSAnimator			animator;
	private int						counter;	// Frame counter

	private final Model				model;

	private final KeyHandler			keyHandler;
	private final MouseHandler			mouseHandler;

	private final Deque<Point2D.Double>			special;
	private final ArrayList<Deque<Point2D.Double>>	regions;

	// special arrays for turning a deque of Point2D.Doubles into a variable size array for use in several methods
	Object arr[] = new Object[12];
	Point2D.Double polyVertices[] = new Point2D.Double[12];

	// Reference Vector
	// TODO: PUT MEMBERS FOR THE REFERENCE VECTOR HERE
	private Point2D.Double refVec;	// reference vector for keeping track of point velocity
	private Point2D.Double p1; // track the points of the first and second vertices of a polygon's side in the for loop
	private Point2D.Double p2;
	private Point2D.Double hp; // the hitpoint of the closest side for the moving point
	private Point2D.Double n; // the normal of the closest side to the moving point

	// Tracer and Bounces
	// TODO: PUT MEMBERS FOR THE TRACER AND BOUNCES HERE
	private ArrayList<Point2D.Double> tracer; // keeps track of object trajectory and renders up to <=1 second
	private ArrayList<Point2D.Double> bounces; // keeps track of object hit points
	private ArrayList<Long> age; // keeps track of system time in milliseconds, for comparing age of elements in tracer


	//**********************************************************************
	// Constructors and Finalizer
	//**********************************************************************

	public View(GLJPanel canvas)
	{
		this.canvas = canvas;

		// Initialize rendering
		counter = 0;
		canvas.addGLEventListener(this);

		// Initialize model (scene data and parameter manager)
		model = new Model(this);

		// Initialize container polygons
		special = createSpecialPolygon();					// For N = 2
		regions = new ArrayList<Deque<Point2D.Double>>();	// For MIN to MAX

		for (int i=MIN_SIDES; i<=MAX_SIDES; i++)
			regions.add(createPolygon(i));

		// Initialize reference vector
		// TODO: INITIALIZE MEMBERS FOR THE REFERENCE VECTOR HERE
		refVec = new Point2D.Double(0.02, -0.01);

		// Initialize tracer and bounces
		// TODO: INITIALIZE MEMBERS FOR THE TRACER AND BOUNCES HERE
		tracer = new ArrayList<Point2D.Double>();
		bounces = new ArrayList<Point2D.Double>();
		age = new ArrayList<Long>();

		// Initialize controller (interaction handlers)
		keyHandler = new KeyHandler(this, model);
		mouseHandler = new MouseHandler(this, model);

		// Initialize animation
		animator = new FPSAnimator(canvas, DEFAULT_FRAMES_PER_SECOND);
		animator.start();
	}

	//**********************************************************************
	// Getters and Setters
	//**********************************************************************

	public GLJPanel	getCanvas()
	{
		return canvas;
	}

	public int	getWidth()
	{
		return w;
	}

	public int	getHeight()
	{
		return h;
	}

	//**********************************************************************
	// Public Methods
	//**********************************************************************

	public void	clearAllTrace()
	{
		// Remove all trajectory and bounce points

		// TODO: YOUR CODE HERE
		age.clear();
		bounces.clear();
		tracer.clear();
	}

	//**********************************************************************
	// Override Methods (GLEventListener)
	//**********************************************************************

	public void	init(GLAutoDrawable drawable)
	{
		w = drawable.getSurfaceWidth();
		h = drawable.getSurfaceHeight();

		renderer = new TextRenderer(new Font("Monospaced", Font.PLAIN, 12),
									true, true);

		initPipeline(drawable);
	}

	public void	dispose(GLAutoDrawable drawable)
	{
		renderer = null;
	}

	public void	display(GLAutoDrawable drawable)
	{
		updatePipeline(drawable);
		update(drawable);
		render(drawable);
	}

	public void	reshape(GLAutoDrawable drawable, int x, int y, int w, int h)
	{
		this.w = w;
		this.h = h;
	}

	//**********************************************************************
	// Private Methods (Rendering)
	//**********************************************************************

	private void	update(GLAutoDrawable drawable)
	{
		counter++;									// Advance animation counter

		Deque<Point2D.Double>	polygon = getCurrentPolygon();
		Point2D.Double			q = model.getObject();

		updatePointWithReflection(polygon, q);
		model.setObjectInSceneCoordinatesAlt(new Point2D.Double(q.x, q.y));

		// Remove old (>1 second) trajectory and bounce points
		boolean oldPurged = false;
		while(!oldPurged) 
		{
			if(age.size() == 0) break;
			if((System.currentTimeMillis()-age.get(0)) > 1000) 
			{
				age.remove(0);
				tracer.remove(0);
			}
			else oldPurged = true;
		}
		/* for(int i = 0; i < age.size(); i++) 
		{
			if((System.currentTimeMillis()-age.get(i)) > 1000) 
			{

			}
		} */

		// TODO: YOUR CODE HERE
	}

	private void	render(GLAutoDrawable drawable)
	{
		GL2	gl = drawable.getGL().getGL2();

		gl.glClear(GL.GL_COLOR_BUFFER_BIT);		// Clear the buffer

		// Draw the scene
		drawMain(gl);								// Draw main content
		drawMode(drawable);						// Draw mode text

		gl.glFlush();								// Finish and display
	}

	//**********************************************************************
	// Private Methods (Pipeline)
	//**********************************************************************

	private void	initPipeline(GLAutoDrawable drawable)
	{
		GL2	gl = drawable.getGL().getGL2();

		gl.glClearColor(0.0f, 0.0f, 0.0f, 0.0f);	// Black background

		// Make points easier to see on Hi-DPI displays
		gl.glEnable(GL2.GL_POINT_SMOOTH);	// Turn on point anti-aliasing
	}

	private void	updatePipeline(GLAutoDrawable drawable)
	{
		GL2			gl = drawable.getGL().getGL2();
		GLU			glu = GLU.createGLU();

		gl.glMatrixMode(GL2.GL_PROJECTION);		// Prepare for matrix xform
		gl.glLoadIdentity();						// Set to identity matrix
		glu.gluOrtho2D(-1.2, 1.2, -1.2, 1.2);		// 2D translate and scale
	}

	//**********************************************************************
	// Private Methods (Scene)
	//**********************************************************************

	private void	drawMode(GLAutoDrawable drawable)
	{
		GL2		gl = drawable.getGL().getGL2();

		renderer.beginRendering(w, h);

		// Draw all text in light gray
		renderer.setColor(0.75f, 0.75f, 0.75f, 1.0f);

		Point2D.Double	cursor = model.getCursor();

		if (cursor != null)
		{
			String		sx = FORMAT.format(new Double(cursor.x));
			String		sy = FORMAT.format(new Double(cursor.y));
			String		s = "Pointer at (" + sx + "," + sy + ")";

			renderer.draw(s, 2, 2);
		}
		else
		{
			renderer.draw("No Pointer", 2, 2);
		}

		String		sn = ("[q|w] Number = " + model.getNumber());
		String		sf = ("[a|s] Factor = " + FORMAT.format(model.getFactor()));
		String		sc = ("[c]   Center moving object in polygon");

		renderer.draw(sn, 2, 16);
		renderer.draw(sf, 2, 30);
		renderer.draw(sc, 2, 44);

		renderer.endRendering();
	}

	private void	drawMain(GL2 gl)
	{
		drawAxes(gl);						// X and Y axes
		drawContainer(gl);					// Container polygon
		drawTracing(gl);					// Object trajectory
		drawBounces(gl);					// Reflection points
		drawObject(gl);					// The moving object
		drawCursor(gl);					// Cursor around the mouse point
	}

	// Draw horizontal (y==0) and vertical (x==0) axes
	private void	drawAxes(GL2 gl)
	{
		gl.glColor3f(0.25f, 0.25f, 0.25f);			// Dark gray

		gl.glBegin(GL.GL_LINES);

		gl.glVertex2d(-10.0, 0.0);
		gl.glVertex2d(10.0, 0.0);

		gl.glVertex2d(0.0, -10.0);
		gl.glVertex2d(0.0, 10.0);

		gl.glEnd();
	}

	// Fills and edges the polygon that is surrounding the moving object.
	private void	drawContainer(GL2 gl)
	{
		Deque<Point2D.Double>	polygon = getCurrentPolygon();

		gl.glColor3f(0.15f, 0.15f, 0.15f);			// Very dark gray
		fillPolygon(gl, polygon);

		gl.glColor3f(1.0f, 1.0f, 1.0f);			// White
		edgePolygon(gl, polygon);
	}

	// If the cursor point is not null, draw something helpful around it.
	private void	drawCursor(GL2 gl)
	{
		Point2D.Double	cursor = model.getCursor();

		if (cursor == null)
			return;

		setColor(gl, 155, 155, 255);		// Bright Yellow
		
		// TODO: YOUR CODE HERE
		gl.glBegin(GL2.GL_LINE_LOOP);

		double	sx = new Double(cursor.x);
		double	sy = new Double(cursor.y);

		gl.glVertex2d(sx, sy-.07);
		gl.glVertex2d(sx+.07, sy);
		gl.glVertex2d(sx, sy+.07);
		gl.glVertex2d(sx-.07, sy);

		gl.glEnd();
	}

	// Draw the moving object, which in this assignment is a single point.
	private void	drawObject(GL2 gl)
	{
		Point2D.Double	object = model.getObject();

		//	TODO: YOUR CODE HERE
		
		/* 
		 * screen borders
		 * 
		 * -x: -1.190
		 * +x:  1.190
		 * -y: -1.190
		 * +y:  1.190
		 */
		
		gl.glPointSize(2.0f);					// Set point size (in pixels)
		setColor(gl, 255, 255, 0);		// Bright Yellow

		// if collision detected w/ polygon, adjust vector trajectory according to reflection angle
		// otherwise, travel normally
		gl.glBegin(GL.GL_POINTS);
		
		gl.glVertex2d(object.x, object.y);

		
		gl.glEnd();
	}

	// Draw the object trajectory in the polygon.
	private void	drawTracing(GL2 gl)
	{
		// TODO: YOUR CODE HERE
		gl.glPointSize(2.0f);				// Set point size (in pixels)
		setColor(gl, 255, 0, 0);		// Red

		// if collision detected w/ polygon, adjust vector trajectory according to reflection angle
		// otherwise, travel normally
		gl.glBegin(GL.GL_POINTS);
		
		for(int i = 0; i < tracer.size(); i++) 
		{
			// System.out.println(tracer.get(0).toString());
			gl.glVertex2d(tracer.get(i).x, tracer.get(i).y);
		}
		
		gl.glEnd();
	}

	// Draw the reflection points on the polygon.
	private void	drawBounces(GL2 gl)
	{
		// TODO: YOUR CODE HERE
		gl.glPointSize(2.0f);				// Set point size (in pixels)
		setColor(gl, 0, 255, 0);		// Green

		// if collision detected w/ polygon, adjust vector trajectory according to reflection angle
		// otherwise, travel normally
		gl.glBegin(GL.GL_POINTS);
		
		for(int i = 0; i < bounces.size(); i++) 
		{
			gl.glVertex2d(bounces.get(i).x, bounces.get(i).y);
		}
		
		gl.glEnd();
	}

	//**********************************************************************
	// Private Methods (Polygons)
	//**********************************************************************

	// Custom polygon for the sides=2 case. Irregular but convex.
	private Deque<Point2D.Double>	createSpecialPolygon()
	{
		Deque<Point2D.Double>	polygon = new ArrayDeque<Point2D.Double>(10);

		polygon.add(new Point2D.Double( 1.00, -0.86));
		polygon.add(new Point2D.Double( 1.00, -0.24));
		polygon.add(new Point2D.Double( 0.48,  0.90));
		polygon.add(new Point2D.Double( 0.05,  1.00));
		polygon.add(new Point2D.Double(-0.34,  0.87));

		polygon.add(new Point2D.Double(-0.86,  0.40));
		polygon.add(new Point2D.Double(-1.00,  0.04));
		polygon.add(new Point2D.Double(-0.93, -0.42));
		polygon.add(new Point2D.Double(-0.53, -0.84));
		polygon.add(new Point2D.Double( 0.71, -1.00));

		// convert polygon points to an object array, save to private variables for use in other methods
		arr = polygon.toArray();
		polyVertices = new Point2D.Double[arr.length];

		for (int i=0; i < arr.length; i++)	// convert object array into a Point2D.Double array
		{
			polyVertices[i] = (Point2D.Double)arr[i];
			// System.out.println(arr2[i].getX() + "\n\n\n\n");
		}

		/* System.out.println("\npX: " + polygon.peek().getX() + "\npY: " + polygon.peek().getY()+ "\n");
		polygon.pop();
		System.out.println("\npX: " + polygon.peek().getX() + "\npY: " + polygon.peek().getY()+ "\n"); */

		// System.out.println("polygon.size(): " + polygon.size());
		return polygon;
	}

	// Creates a regular N-gon with points stored in counterclockwise order.
	// The polygon is centered at the origin with first vertex at (1.0, 0.0).
	public Deque<Point2D.Double>	createPolygon(int sides)
	{
		Deque<Point2D.Double>	polygon = new ArrayDeque<Point2D.Double>(sides);
		double x = 1.00;
		double y = 1.00;
		double polyTheta = 0.0 * Math.PI;

		// TODO: YOUR CODE HERE
		for(int i = 0; i < sides; i++)
		{
			double polyX = x*Math.cos(polyTheta);
			double polyY = y*Math.sin(polyTheta);

			polygon.add(new Point2D.Double( polyX, 
			polyY));
			polyTheta = polyTheta + 2.0*Math.PI/sides;
		}
		return polygon;
	}

	// Draws the sides of the specified polygon.
	private void	edgePolygon(GL2 gl, Deque<Point2D.Double> polygon)
	{
		// TODO: YOUR CODE HERE
		
		gl.glBegin(GL2.GL_LINE_LOOP);

		for (int i=0; i < arr.length; i++)
		{
			polyVertices[i] = (Point2D.Double)arr[i];
			// System.out.println(arr2[i].getX() + "\n\n\n\n");
		}
		for (int i=0; i < polyVertices.length; i++)
		{
			gl.glVertex2d(polyVertices[i].getX(), polyVertices[i].getY());
		}

		gl.glEnd();
	}

	// Draws the interior of the specified polygon.
	private void	fillPolygon(GL2 gl, Deque<Point2D.Double> polygon)
	{
		// TODO: YOUR CODE HERE
		
		gl.glBegin(GL2.GL_POLYGON);

		// convert polygon points to an object array, save to private variables for use in other methods
		arr = polygon.toArray();
		polyVertices = new Point2D.Double[arr.length];

		for (int i=0; i < arr.length; i++)	// convert object array into a Point2D.Double array
		{
			polyVertices[i] = (Point2D.Double)arr[i];
			// System.out.println(arr2[i].getX() + "\n\n\n\n");
		}
		for (int i=0; i < polyVertices.length; i++)	 // draw each point until looping back around
		{
			gl.glVertex2d(polyVertices[i].getX(), polyVertices[i].getY());
		}

		gl.glEnd();
	}

	// Get the polygon that is currently containing the moving object.
	private Deque<Point2D.Double>	getCurrentPolygon()
	{
		int	sides = model.getNumber();

		if (sides == 2)
		{
			return special;
		}
		else if ((MIN_SIDES <= sides) && (sides <= MAX_SIDES))
		{
			return regions.get(sides - MIN_SIDES);
		}
		else
			return null;
	}

	// Special method for privileged use by the Model class ONLY.
	public boolean	currentPolygonContains(Point2D.Double q)
	{
		return contains(q);
	}

	//**********************************************************************
	// Private Methods (Reflection)
	//**********************************************************************

	// Updates the x and y coordinates of point q. Adds a vector to the provided
	// point, reflecting as needed off the sides of the provided polygon to
	// determine the new coordinates. The new coordinates are "returned" in q.
	public void	updatePointWithReflection(Deque<Point2D.Double> polygon,
											  Point2D.Double q)
	{
		// TODO: YOUR CODE HERE. Hints for how to approach it follow.
		Point2D.Double object = model.getObject();

		// Use the reference vector to remember the current direction of
		// movement with a magnitude equal to the default distance (factor=1.0).
		// For each update, copy the reference vector and scale it by the
		// current speed factor...
		Point2D.Double refScale = new Point2D.Double(refVec.x*model.getFactor(), refVec.y*model.getFactor());

		// 1. Calculate which side the point will reach first.

		// variables for interpolating the side with a point on it that the moving point will reach first
		hp = null;
		n = null;
		double sd = 99999999; // narrows down the side with the closest distance to the moving point

		// Loop the polygon counterclockwise, taking vertices pairwise.
		// For each side, see "Intersection of a Line through a Line".
		// Important: Check for edge cases (pun?) in which q is parallel
		// to the side or slightly outside it (due to roundoff error).
		// See Figure 4.37 and the dot products below it on page 176.
		// Always remember to check for divide-by-zero!
		// 
		// P_hit = R + v((n dot (Q-R))/(n dot v))
		// P = R + vt		<-- R is starting point, aka the 'b' in y = mx + b
		// v is reference vector, unaltered by point speed which goes from 0.0 to 1.0
		// in context of view.object, this makes R the location of object when no point speed is applied
		// n is (-(p2-p1).y, (p2-p1).x)
		// Q is the first of the two vertices of any side, i.e. the q_j followed by q_j+1
		//
		// to figure out which side the point will hit first before calulating the trajectory, all we need to do is the following steps:
		// 1. find dot product of vector and each side's normal (-y, x)
		// 2. the side the vector will reach first *should* be the side with the smallest positive dot product
		// 3. plug said side into reflecting trajectories bottom left equations and go from there

		// loop through each vertex and find the one that produces the side the moving point will hit next
		for(int i = 0; i < polyVertices.length; i++)
		{
			if(i == polyVertices.length-1) // if last vertex, pair with 0th index
			{
				// pair vertex on current index with that of following index
				p1 = new Point2D.Double(polyVertices[i].x, polyVertices[i].y);
				p2 = new Point2D.Double(polyVertices[0].x, polyVertices[0].y);
			}
			else 
			{
				p1 = new Point2D.Double(polyVertices[i].x, polyVertices[i].y);
				p2 = new Point2D.Double(polyVertices[i+1].x, polyVertices[i+1].y);
			}

			Point2D.Double vo = new Point2D.Double(object.x-p1.x, object.y-p1.y); // find vector between first vertex and object...
			Point2D.Double vi = new Point2D.Double(p2.x-p1.x, p2.y-p1.y); // ...then find vector from vertices...
			//Point2D.Double ni = new Point2D.Double(-p.y, p.x); // ...then find normal from vector

			// projection a = v . n / n . n
			double a = dot(vo.x, vo.y, 0, vi.x, vi.y, 0)/dot(vi.x, vi.y, 0, vi.x, vi.y, 0);

			// make the projection a unit
			if (a>1) a = 1;
			else if (a<0) a = 0;

			// use interpolation to find point along side closest to the moving object
			Point2D.Double ip = new Point2D.Double(p1.x + a*vi.x, p1.y + a*vi.y);

			// find distance between object and interpolated point
			double d = ip.distance(object);

			if(d < sd) // becomes new closest side if smaller distance than current smallest side
			{
				double l = Math.sqrt(dot(vi.x, vi.y, 0, vi.x, vi.y, 0)); // length of side
				n = new Point2D.Double(-vi.y/l, vi.x/l); // normal of new closest side, scaled to length
				hp = ip; // set hitpoint of object on side to the interpolated point
				sd = d; // finally, update closest distance to object
			}
		}

		// 2. If the point WON'T reach the closest side in this update,
		//    simply add the vector to it, and break out of the loop.
		//    Or if it WILL reach the side:
		//    Move the point to the hit point on the closest side.
		//    Calculate the reflected vector. Reduce it by the amount
		//    already moved, and use it as the new scaled vector.
		//    After bounces, remember to update the direction of the reference
		//    vector for the next time updatePointWithReflection() is called!

		// the object is contained after it is moved, then move it no problem
		if(contains(new Point2D.Double(object.x+refScale.x, object.y+refScale.y)) == true) 
		{
			// add new point to object trail and keep track of its add time
			age.add(System.currentTimeMillis());
			tracer.add(new Point2D.Double(object.x, object.y));

			// update position of object
			model.setObjectInSceneCoordinates(new Point2D.Double(object.x+refScale.x, object.y+refScale.y));
		}
		else // otherwise, implement reflection
		{
			// add new bounce point to list
			bounces.add(new Point2D.Double(hp.x, hp.y));

			// move object to interpolated hit point
			model.setObjectInSceneCoordinates(new Point2D.Double(hp.x, hp.y));

			// update scaled vector
			refScale.x -= 2 * dot(refScale.x, refScale.y, 0, n.x, n.y, 0) * n.x;
			refScale.y -= 2 * dot(refScale.x, refScale.y, 0, n.x, n.y, 0) * n.y;

			// determine ratio of original speed to new speed
			double rx = Math.abs(refVec.x);
			double ry = Math.abs(refVec.y);
			double rsx = Math.abs(refScale.x/model.getFactor());
			double rsy = Math.abs(refScale.y/model.getFactor());
			double ratio = Math.abs((rx+ry)/(rsx+rsy));

			/* System.out.println("\nrefVec: " + refVec.toString()
								+ "\nrefScale: " + refScale.toString()); */

			// change current vector to new direction, maintain speed between vectors
			refVec.x = (refScale.x/model.getFactor())*ratio;
			refVec.y = (refScale.y/model.getFactor())*ratio;
		}
	}

	//**********************************************************************
	// Private Methods (Vectors)
	//**********************************************************************

	// This might be a method to calculate a dot product. Sure seems like it.
	private double		dot(double vx, double vy, double vz,
							double wx, double wy, double wz)
	{
		// TODO: YOUR CODE HERE

		double dot = vx*wx + vy*wy + vz*wz;
		return dot;
	}

	// Determines if point q is to the left of line p1->p2. If strict is false,
	// points exactly on the line are considered to be left of it.
	private boolean	isLeft(Point2D.Double p1, Point2D.Double p2,
							   Point2D.Double q, boolean strict)
	{
		// TODO: YOUR CODE HERE

		// Hint: Use dot(). See the slide on "Testing Containment in 2D".

		Point2D.Double sp = new Point2D.Double(p2.x - p1.x, p2.y - p1.y); // finding side vector

		Point2D.Double dif_q = new Point2D.Double(q.x - p1.x, q.y - p1.y); // distance from point to first vertex of side
		Point2D.Double dif_p = new Point2D.Double(-sp.y, sp.x);	// normal side vector
		double dp = dot(dif_q.x, dif_q.y, 0, dif_p.x, dif_p.y, 0); // dot of the two

		if(dp > 0 || (dp == 0 && !strict)) 
			return true;
		else 
			return false;
	}

	// Determines if point q is inside a polygon. The polygon must be convex
	// with points stored in counterclockwise order. Points exactly on any side
	// of the polygon are considered to be outside of it.
	private boolean	contains(Point2D.Double q)
	{
		// TODO: YOUR CODE HERE

		// Hint: Use isLeft(). See the slide on "Testing Containment in 2D".

		for(int i = 0; i < polyVertices.length; i++) 
		{
			if(i == polyVertices.length-1) 
			{
				if(!isLeft(polyVertices[i], polyVertices[0], q, true)) 
				{
					return false;
				}
			}
			else 
			{
				if(!isLeft(polyVertices[i], polyVertices[i+1], q, true)) 
				{
					return false;
				}
			}
		}
		return true;
	}


// utility functions for setting color
private void	setColor(GL2 gl, int r, int g, int b, int a)
	{
		gl.glColor4f(r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f);
	}

	private void	setColor(GL2 gl, int r, int g, int b)
	{
		setColor(gl, r, g, b, 255);
	}
}
//******************************************************************************
