#include "DrawFluid.h"
#include "glad/glad.h"

DrawFluid::DrawFluid(std::shared_ptr<Fluid> inFluid)
{
	FromFluid(inFluid);
}

void DrawFluid::FromFluid(std::shared_ptr<Fluid> inFluid)
{
	mParticlePoints.reserve(inFluid.get()->GetNumParticles());

	std::vector<Particle>::iterator itr;

	for (itr = inFluid.get()->mParticles.begin(); itr < inFluid.get()->mParticles.end(); itr++)
	{
		mParticlePoints.push_back(itr->GetPosition());
	}
}

void DrawFluid::UpdateOpenGLBuffers()
{

}

