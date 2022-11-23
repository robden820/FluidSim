#include "DrawFluid.h"
#include "glad/glad.h"

DrawFluid::DrawFluid(std::shared_ptr<Fluid> inFluid)
{
	FromFluid(inFluid);
}

void DrawFluid::FromFluid(std::shared_ptr<Fluid> inFluid)
{
	mParticlePoints.clear();
	mParticlePoints.reserve(inFluid->GetNumParticles());

	std::vector<Particle>::iterator itr;

	for (itr = inFluid->GetParticles().begin(); itr < inFluid->GetParticles().end(); itr++)
	{
		mParticlePoints.push_back(itr->GetPosition());
	}
}

void DrawFluid::UpdateOpenGLBuffers()
{

}

