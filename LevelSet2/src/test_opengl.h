int width , height;
const GLfloat vertex[] = {
	-0.9 , -0.9 , -2 , 0.9 , -0.9 , -2 , 0 , 0.9 , -2
};

void Draw(void);
void Draw() {
	glDrawArrays(GL_POLYGON , 0 , 3);
	glViewport(0 , 0 , width , height);
}

void disp( void ) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
	glEnable(GL_DEPTH_TEST);

	glPushMatrix();
		glColor3f(1 , 0 ,0);
		glTranslatef(-0.5 , 0 , 0);
		Draw();
	glPopMatrix();
	glPushMatrix();
		glColor3f(0 , 0 , 1);
		glTranslatef(0.5 , 0 , -1);
		Draw();
	glPopMatrix();

	glFlush();
}

void reshape(int w , int h) {
	width = w; height = h;
	disp();
}

int test_opengl(int argc , char ** argv) {
	glutInit(&argc , argv);
	glutInitWindowPosition(100 , 50);
	glutInitWindowSize(400 , 300);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH);

	glutCreateWindow("Kitty on your lap");
	glutDisplayFunc(disp);
	glutReshapeFunc(reshape);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-1 , 1 , -1 , 1 , 2 , 10);

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3 , GL_FLOAT , 0 , vertex);

	glutMainLoop();
	return 0;
}
